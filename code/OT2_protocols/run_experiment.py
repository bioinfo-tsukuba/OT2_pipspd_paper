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
        for i in range(mix_before[0]):
            pipette.mix(1, mix_before[1], source)

    # transfering
    pipetting_record.append_source(
        pipette,
        source,
        volume,
        remark
    )
    pipetting_record.append_dest(
        pipette,
        dest,
        volume,
        remark
    )
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
        for i in range(mix_after[0]):
            pipette.mix(1, mix_after[1], dest)

    if drop_tip:
        pipette.drop_tip()


# This function is for calculation of agar height of microplate 8x12 #################
def agar_height_from_weight(source_csv, deck_num):
    if (deck_num==6):
        value_loc = 0
    elif (deck_num==9):
        value_loc = 1
    # DISH_THICKNESS = 3  # mm
    DISH_THICKNESS = 2  # mm
    DISH_WEIGHT = 35.28
    AGAR_DENSITY = 0.00102  # g/mm^3
    DISH_BOTTOM_AREA = 9458.02 # mm^2
    agar_plate_weight = pd.read_csv(source_csv, header=None)
    agar_plate_weight = agar_plate_weight.values[value_loc][1]
    agar_weight = agar_plate_weight - DISH_WEIGHT
    agar_height = round(
        agar_weight / (DISH_BOTTOM_AREA * AGAR_DENSITY), 1
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
    'protocolName': 'Spotting test for microplate 8x12',
    'description': '''Nothing''',
    'author': 'Shodai Taguchi'
}


def run(protocol: protocol_api.ProtocolContext):

    # hyperparameters #########################################################

    USE_SINGLE_TECHREP = False
    COLOR1_POSITION = 'A3'
    COLOR2_POSITION = 'A4'
    COLOR3_POSITION = 'B4'
    WATER_POSITION = 'B3'
    TIP300_FIRST_TIP = 'E2'
    TIP20_FIRST_TIP = 'A1'
    # MIX_VOLUME = 200
    # MIX_TIMES = 3


    # Config about spotting
    SPOT_SPEED = 50
    SPOT_VOLUME = 2
    # SPOT_HEIGHT = agar_height_from_weight('./agar_plate_weight.csv')
    SPOT_HEIGHT_6 = agar_height_from_weight('./agar_plate_weight.csv', 6)
    SPOT_HEIGHT_9 = agar_height_from_weight('./agar_plate_weight.csv', 9)

    # labware definition ######################################################

    # Load original labware definition
    def_microplate = load_definition_from_json(
        "./labware_OT2_yeast/"
        + "tgc_spotassay_microplate_8x12.json"
    )
    
    # load labware ############################################################

    # Tip rack
    tip300 = protocol.load_labware('opentrons_96_tiprack_300ul', 4)
    tip20 = protocol.load_labware('opentrons_96_tiprack_20ul', 5)

    # Tube rack
    # Source of Colored Water
    tube_source = protocol.load_labware(
        'opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical', 7
    )

    # 96 well plate for transport system gen2
    plate_96well = protocol.load_labware('corning_96_wellplate_360ul_flat', 8)

    # Agar plate for spot-assay
    # agar_plate = protocol.load_labware_from_definition(def_petri_dish, 9)
    agar_plate_6 = protocol.load_labware_from_definition(def_microplate, 6)
    agar_plate_9 = protocol.load_labware_from_definition(def_microplate, 9)



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
    plate_96well.set_offset(-3.20, 0.60, 0.00)
    p300.starting_tip = tip300.well(TIP300_FIRST_TIP)
    p20.starting_tip = tip20.well(TIP20_FIRST_TIP)

    pipetting_record = PipettingRecordList()

    # Alias to each source
    source_color1 = tube_source[COLOR1_POSITION]
    source_color2 = tube_source[COLOR2_POSITION]
    source_color3 = tube_source[COLOR3_POSITION]
    water = tube_source[WATER_POSITION]

    # =========================================================================

    # [index, biorep, speed_index, remark]
    dilution_order = [
        (0,  source_biorep1, 3, "B1_S4_T1"),
        (1,  source_biorep1, 2, "B1_S3_T1"),
        (2,  source_biorep1, 1, "B1_S2_T1"),
        (3,  source_biorep1, 0, "B1_S1_T1"),
        (4,  source_biorep2, 3, "B2_S4_T1"),
        (5,  source_biorep2, 2, "B2_S3_T1"),
        (6,  source_biorep2, 1, "B2_S2_T1"),
        (7,  source_biorep2, 0, "B2_S1_T1"),
        (8,  source_biorep1, 3, "B1_S4_T2"),
        (9,  source_biorep1, 2, "B1_S3_T2"),
        (10, source_biorep1, 1, "B1_S2_T2"),
        (11, source_biorep1, 0, "B1_S1_T2"),
        (12, source_biorep2, 3, "B2_S4_T2"),
        (13, source_biorep2, 2, "B2_S3_T2"),
        (14, source_biorep2, 1, "B2_S2_T2"),
        (15, source_biorep2, 0, "B2_S1_T2"),
        (16, source_biorep1, 3, "B1_S4_T3"),
        (17, source_biorep1, 2, "B1_S3_T3"),
        (18, source_biorep1, 1, "B1_S2_T3"),
        (19, source_biorep1, 0, "B1_S1_T3"),
        (20, source_biorep2, 3, "B2_S4_T3"),
        (21, source_biorep2, 2, "B2_S3_T3"),
        (22, source_biorep2, 1, "B2_S2_T3"),
        (23, source_biorep2, 0, "B2_S1_T3"),
    ]

    # Before yeast transfer, add water to each well

    # Add water to RNA-seq tube
    aspirated_color1 = 0
    aspirated_color2 = 0
    aspirated_color3 = 0
    aspirated_water = 0
    for i in range(96):
        surface_height_decrese_color1 = volume_to_height_50mL_tube(100)
        surface_height_decrese_color2 = volume_to_height_50mL_tube(100)
        surface_height_decrese_color3 = volume_to_height_50mL_tube(100)
        surface_height_decrese_water = volume_to_height_50mL_tube(100)
        
        if i % 3 == 0:
            transfer_with_timestamp(
                p300,
                100,
                source_color1.bottom(z=95 - aspirated_color1),
                plate_96well.wells()[i],
                pipetting_record,
                drop_tip=True,
                remark="Add color1 to 96 well plate"
            )
            aspirated_color1 += surface_height_decrese_color1
        
        elif i % 3 == 1:
            transfer_with_timestamp(
                p300,
                100,
                source_color2.bottom(z=95 - aspirated_color2),
                plate_96well.wells()[i],
                pipetting_record,
                drop_tip=True,
                remark="Add color2 to 96 well plate"
            )
            aspirated_color2 += surface_height_decrese_color2
        elif i % 3 == 2:
            transfer_with_timestamp(
                p300,
                100,
                source_color3.bottom(z=95 - aspirated_color3),
                plate_96well.wells()[i],
                pipetting_record,
                drop_tip=True,
                remark="Add color3 to 96 well plate"
            )
            aspirated_color3 += surface_height_decrese_color3
    
    agars = [agar_plate_6, agar_plate_9]
    heights = [SPOT_HEIGHT_6, SPOT_HEIGHT_9]
    # Spotting on Agar in 96 well plate
    for j in range(96):
        for agar_plate_label, agar_plate_height in agars, heights:
            agar_plate_row_flatten = list_2d_flatten(agar_plate_label.rows())
            
            if j % 3 == 0:
                transfer_with_timestamp(
                    p20,
                    SPOT_VOLUME,
                    plate_96well.wells()[j],
                    agar_plate_row_flatten[j].bottom(z=agar_plate_height),
                    pipetting_record=pipetting_record,
                    drop_tip=False,
                    remark=agar_plate_label+"Color1_spotting"
                    )
            
            elif j % 3 == 1:
                transfer_with_timestamp(
                    p20,
                    SPOT_VOLUME,
                    plate_96well.wells()[j],
                    agar_plate_row_flatten[j].bottom(z=SPOT_HEIGHT_6),
                    pipetting_record=pipetting_record,
                    drop_tip=False,
                    remark=agar_plate_label+"Color2_spotting"
                )
            
            elif j % 3 == 2:
                transfer_with_timestamp(
                    p20,
                    SPOT_VOLUME,
                    plate_96well.wells()[j],
                    agar_plate_row_flatten[j].bottom(z=SPOT_HEIGHT_6),
                    pipetting_record=pipetting_record,
                    drop_tip=False,
                    remark=agar_plate_label+"Color3_spotting"
                    )
        p20.drop_tip()

    # Output record
    pipetting_record.export_csv('pipetting_record.csv')
