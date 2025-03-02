from opentrons import protocol_api
from opentrons.types import Location
import json
import time
import pandas as pd
import math
import itertools


class PipettingRecord():
    def __init__(
            self,
            pipette: protocol_api.InstrumentContext,
            well,
            operation: str,
            volume: float,
            remark: str = ''
    ) -> None:
        if type(well) == Location:
            self.well = str(well.labware)
        else:
            self.well = str(well)
        self.operation = operation
        self.volume = str(volume)
        self.unix_time = time.time()
        self.remark = remark
        self.tip_from = pipette._last_tip_picked_up_from

    def __str__(self) -> str:
        return (
            '"' + self.well + '"'
            + ','
            + self.operation
            + ','
            + self.volume
            + ','
            + str(self.unix_time)
            + ','
            + str(to_JST(self.unix_time))
            + ','
            + str(self.tip_from)
            + ','
            + self.remark
        )


# Load labware definition from JSON file
def load_definition_from_json(file_path):
    with open(file_path) as file:
        definition = json.load(file)
    return definition


# Convert unix time to JST (GMT+9) regardless of the system timezone
def to_JST(unix_time):
    JST_OFFSET = 9 * 60 * 60
    jst_time = time.gmtime(unix_time + JST_OFFSET)
    return time.strftime('%Y-%m-%d %H:%M:%S', jst_time)


def change_pipette_speed(
        pipette: protocol_api.InstrumentContext,
        aspirate_speed,
        dispense_speed
):
    pipette.flow_rate.aspirate = aspirate_speed
    pipette.flow_rate.dispense = dispense_speed


class PipettingRecordList():
    def __init__(self) -> None:
        self.pipetting_record = []

    def append_mix(self, pipette, well, volume, remark=''):
        self.pipetting_record.append(
            PipettingRecord(pipette, well, 'mix', volume, remark)
        )

    def append_source(self, pipette, well, volume, remark=''):
        self.pipetting_record.append(
            PipettingRecord(pipette, well, 'source', volume, remark)
        )

    def append_dest(self, pipette, well, volume, remark=''):
        self.pipetting_record.append(
            PipettingRecord(pipette, well, 'dest', volume, remark)
        )

    def export_csv(self, file_path):
        pipetting_record_file = open(file_path, 'wt')
        print(
            'well,operation,volume,unix_time,JST(GMT+9),tip_from,remark',
            file=pipetting_record_file
        )
        [
            print(
                self.pipetting_record[i],
                file=pipetting_record_file
            ) for i in range(len(self.pipetting_record))
        ]
        pipetting_record_file.close


def transfer_with_timestamp(
        pipette: protocol_api.InstrumentContext,
        volume: float,
        source: protocol_api.labware.Well,
        dest: protocol_api.labware.Well,
        pipetting_record: PipettingRecordList,
        remark: str = '',
        mix_before: tuple = None,
        mix_after: tuple = None,
        drop_tip: bool = True
):
    # check if pipette has tip
    if not pipette.has_tip:
        pipette.pick_up_tip()

    # pre-transfering mixing
    if not (mix_before is None):
        # mix_before is specified
        pipetting_record.append_mix(
            pipette,
            source,
            mix_before[1],
            remark
        )
        # pipetting_record.append(
        #     PipettingRecord(pipette, source, 'mix', mix_before[1], remark)
        # )
        for i in range(mix_before[0]):
            pipette.mix(1, mix_before[1], source)

    # transfering
    pipetting_record.append_source(
        pipette,
        source,
        volume,
        remark
    )
    # pipetting_record.append(
    #     PipettingRecord(pipette, source, 'source', volume, remark)
    # )
    pipetting_record.append_dest(
        pipette,
        dest,
        volume,
        remark
    )

    # pipetting_record.append(
    #     PipettingRecord(pipette, dest, 'dest', volume, remark)
    # )
    pipette.transfer(volume, source, dest, new_tip='never')

    # post-transfering mixing
    if not (mix_after is None):
        # mix_after is specified
        pipetting_record.append_mix(
            pipette,
            dest,
            mix_after[1],
            remark
        )
        # pipetting_record.append(
        #     PipettingRecord(pipette, dest, 'mix', mix_after[1], remark)
        # )
        for i in range(mix_after[0]):
            pipette.mix(1, mix_after[1], dest)

    if drop_tip:
        pipette.drop_tip()


def agar_height_from_weight(source_csv, deck_num):
    source_csv = pd.read_csv(source_csv, header=None)
    if (deck_num==6):
        value_loc = 0
    elif (deck_num==9):
        value_loc = 1
    
    plate_type = source_csv.values[value_loc][2]
    
    if (plate_type == "microplate"):
        DISH_WEIGHT = 35.28 # g
        AGAR_DENSITY = 0.00102  # g/mm^3
        PLATE_BOTTOM_AREA = 9856.785 # mm^2
        
        agar_plate_weight = source_csv.values[value_loc][1]    
        agar_weight = agar_plate_weight - DISH_WEIGHT
        agar_height = round(
            agar_weight / (PLATE_BOTTOM_AREA * AGAR_DENSITY), 1
        )  # calculate native agar height(cm)
        return agar_height
        # plus the bottom thickness of dish (3.80 mm) and make array data for adapting some different dishes
        
    elif (plate_type == "petridish"):
        DISH_THICKNESS = 2  # mm
        DISH_WEIGHT = 17.88
        AGAR_DENSITY = 0.00102  # g/mm^3
        DISH_DIAMETER = 86
        agar_plate_weight = source_csv.values[value_loc][1]
        agar_weight = agar_plate_weight - DISH_WEIGHT
        agar_height = round(
            agar_weight / ((math.pi * (DISH_DIAMETER / 2) ** 2) * AGAR_DENSITY), 1
        )  # calculate native agar height(cm)
        return agar_height + DISH_THICKNESS
        # agar_height_array.append(agar_height+3)  #plus the bottom thickness of
        # dish (3mm) and make array data for adapting some different dishes


def volume_to_height_50mL_tube(water_volume):
    water_height = water_volume / 500
    return water_height


# CAUTION!
# This function only supports 2D array with all the same type of elements.
# e.g. [[1,2,3],[4,5,6]] OK
#      [[1,2,3],[4,5]] NG
def list_2d_flatten(target_list):
    return list(itertools.chain.from_iterable(target_list))


metadata = {
    'apiLevel': '2.13',
    'protocolName': 'yeast RNA-seq (first-half)',
    'description': '''RNA-seq of yeast in iLab''',
    'author': 'Ryosuke Matsuzawa'
}


def run(protocol: protocol_api.ProtocolContext):

    # hyperparameters #########################################################

    USE_SINGLE_TECHREP = False
    source_biorep1_POSITION = 'A4'
    source_biorep2_POSITION = 'B4'
    WATER_POSITION = 'B3'
    TIP300_FIRST_TIP = 'A1'
    TIP20_FIRST_TIP = 'A1'
    # flow rate: [speed1, speed2, speed3, speed4]
    PIPETTE_FLOW_RATE = [50, 130, 210, 290]
    # PIPETTE_FLOW_RATE = [290, 290, 290, 290]
    MIX_VOLUME = 200
    MIX_TIMES = 3

    # labware definition ######################################################

    def_96well_plate = load_definition_from_json(
        "./labware_OT2_yeast/"
        + "corning_96_wellplate_360ul_flat_transportation_unit_ver2.json"
    )

    def_microplate = load_definition_from_json(
        "./labware_OT2_yeast/"
        + "tgc_spotassay_microplate_8x12.json"
    )

    # load labware ############################################################

    # Tip rack
    tip300 = protocol.load_labware('opentrons_96_tiprack_300ul', 4)
    tip20 = protocol.load_labware('opentrons_96_tiprack_20ul', 5)

    # Tube rack
    # Source of yeast & water.
    tube_source = protocol.load_labware(
        'opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical', 7
    )

    # RNA-seq sample & source of spotting
    # It is not used in first-half, but let OT-2 know that it exists.
    tube_RNAseq = protocol.load_labware(
        'opentrons_24_tuberack_nest_2ml_screwcap', 8
    )

    # Agar plate for spot-assay
    # agar_plate = protocol.load_labware_from_definition(def_petri_dish, 9)
    ## agar_plate_6 = protocol.load_labware_from_definition(def_petri_dish, 6)
    ## agar_plate_9 = protocol.load_labware_from_definition(def_petri_dish, 9)
    agar_plate_6 = protocol.load_labware_from_definition(def_microplate, 6)
    agar_plate_9 = protocol.load_labware_from_definition(def_microplate, 9)

    # 96 well plate for transport system gen2
    plate_96well = protocol.load_labware_from_definition(
        def_96well_plate,
        11
    )

    # pipete
    p300 = protocol.load_instrument(
        'p300_single_gen2',
        'left',
        tip_racks=[tip300]
    )
    p20 = protocol.load_instrument(
        'p20_single_gen2',
        'right',
        tip_racks=[tip20]
    )

    # labware config ##########################################################

    # Labware offset
    tip300.set_offset(0.00, 0.70, 0.00)
    plate_96well.set_offset(-1.20, 0.60, 0.00)

    # starting tip
    p300.starting_tip = tip300.well(TIP300_FIRST_TIP)
    p20.starting_tip = tip20.well(TIP20_FIRST_TIP)

    # Alias to each source
    source_biorep1 = tube_source[source_biorep1_POSITION]
    source_biorep2 = tube_source[source_biorep2_POSITION]
    water = tube_source[WATER_POSITION]

    pipetting_record = PipettingRecordList()

    # =========================================================================

    # recipe for absorbance measurement
    # recipe: (source, speed)
    recipe_single_techrep = [
        (source_biorep1, PIPETTE_FLOW_RATE[3]),  # 0
        (source_biorep1, PIPETTE_FLOW_RATE[2]),  # 2
        (source_biorep1, PIPETTE_FLOW_RATE[1]),  # 4
        (source_biorep1, PIPETTE_FLOW_RATE[0]),  # 6

        (source_biorep2, PIPETTE_FLOW_RATE[3]),  # 1
        (source_biorep2, PIPETTE_FLOW_RATE[2]),  # 3
        (source_biorep2, PIPETTE_FLOW_RATE[1]),  # 5
        (source_biorep2, PIPETTE_FLOW_RATE[0]),  # 7
    ]

    recipe_triple_techrep = [
        # TechRep1
        (source_biorep1, PIPETTE_FLOW_RATE[3]),  # 0
        (source_biorep1, PIPETTE_FLOW_RATE[2]),  # 2
        (source_biorep1, PIPETTE_FLOW_RATE[1]),  # 4
        (source_biorep1, PIPETTE_FLOW_RATE[0]),  # 6

        (source_biorep2, PIPETTE_FLOW_RATE[3]),  # 1
        (source_biorep2, PIPETTE_FLOW_RATE[2]),  # 3
        (source_biorep2, PIPETTE_FLOW_RATE[1]),  # 5
        (source_biorep2, PIPETTE_FLOW_RATE[0]),  # 7

        # TechRep2
        (source_biorep1, PIPETTE_FLOW_RATE[3]),  # 0
        (source_biorep1, PIPETTE_FLOW_RATE[2]),  # 2
        (source_biorep1, PIPETTE_FLOW_RATE[1]),  # 4
        (source_biorep1, PIPETTE_FLOW_RATE[0]),  # 6

        (source_biorep2, PIPETTE_FLOW_RATE[3]),  # 1
        (source_biorep2, PIPETTE_FLOW_RATE[2]),  # 3
        (source_biorep2, PIPETTE_FLOW_RATE[1]),  # 5
        (source_biorep2, PIPETTE_FLOW_RATE[0]),  # 7

        # TechRep3
        (source_biorep1, PIPETTE_FLOW_RATE[3]),  # 0
        (source_biorep1, PIPETTE_FLOW_RATE[2]),  # 2
        (source_biorep1, PIPETTE_FLOW_RATE[1]),  # 4
        (source_biorep1, PIPETTE_FLOW_RATE[0]),  # 6

        (source_biorep2, PIPETTE_FLOW_RATE[3]),  # 1
        (source_biorep2, PIPETTE_FLOW_RATE[2]),  # 3
        (source_biorep2, PIPETTE_FLOW_RATE[1]),  # 5
        (source_biorep2, PIPETTE_FLOW_RATE[0]),  # 7
    ]

    # protocol ###############################################################

    well_number_96well = 0
    if USE_SINGLE_TECHREP:
        recipe = recipe_triple_techrep
    else:
        recipe = recipe_single_techrep

    for source, speed in recipe:
        change_pipette_speed(p300, speed, speed)
        transfer_with_timestamp(
            p300,
            200,
            source,
            plate_96well.wells()[well_number_96well],
            pipetting_record,
            mix_before=(MIX_TIMES, MIX_VOLUME)
        )
        well_number_96well += 1

    # =========================================================================
    pipetting_record.export_csv('pipetting_record_first.csv')
