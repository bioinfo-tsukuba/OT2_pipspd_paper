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
