from opentrons import protocol_api
from opentrons.types import Location
import json
import time
import pandas as pd
import math
import itertools


class PipettingRecord():
    class PipettingRecord:
        """
        A class to record pipetting operations performed by a pipette.
        Attributes:
            well (str): The well or labware location where the operation is performed.
            operation (str): The type of operation performed (e.g., "aspirate", "dispense").
            volume (str): The volume of liquid involved in the operation.
            unix_time (float): The Unix timestamp when the operation was recorded.
            remark (str): Additional remarks or comments about the operation.
            tip_from (str): The location from which the pipette tip was last picked up.
        Methods:
            __str__(): Returns a string representation of the pipetting record, including
                       well, operation, volume, timestamps, tip location, and remarks.
        """
    
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


def load_definition_from_json(file_path):
    """
    Load a JSON file and return its contents as a Python dictionary.
    Args:
        file_path (str): The path to the JSON file to be loaded.
    Returns:
        dict: The contents of the JSON file as a dictionary.
    Raises:
        FileNotFoundError: If the specified file does not exist.
        json.JSONDecodeError: If the file is not a valid JSON.
    """

    with open(file_path) as file:
        definition = json.load(file)
    return definition


def to_JST(unix_time):
    """
    Converts a Unix timestamp to a formatted date and time string in Japan Standard Time (JST).
    Args:
        unix_time (int): The Unix timestamp (seconds since epoch) to be converted.
    Returns:
        str: The formatted date and time string in JST (YYYY-MM-DD HH:MM:SS).
    """
    
    JST_OFFSET = 9 * 60 * 60
    jst_time = time.gmtime(unix_time + JST_OFFSET)
    return time.strftime('%Y-%m-%d %H:%M:%S', jst_time)


def change_pipette_speed(
        pipette: protocol_api.InstrumentContext,
        aspirate_speed,
        dispense_speed
):
    """
    Adjusts the flow rates for aspirating and dispensing on a given pipette.
    Parameters:
        pipette (protocol_api.InstrumentContext): The pipette object whose flow rates are to be modified.
        aspirate_speed (float): The desired flow rate for aspirating (in µL/s).
        dispense_speed (float): The desired flow rate for dispensing (in µL/s).
    Returns:
        None
    """
    pipette.flow_rate.aspirate = aspirate_speed
    pipette.flow_rate.dispense = dispense_speed


class PipettingRecordList():
    """
    A class to manage and record pipetting operations.
    Methods
    -------
    __init__():
        Initializes an empty list to store pipetting records.
    append_mix(pipette, well, volume, remark=''):
        Appends a mixing operation record to the pipetting record list.
    append_source(pipette, well, volume, remark=''):
        Appends a source operation record to the pipetting record list.
    append_dest(pipette, well, volume, remark=''):
        Appends a destination operation record to the pipetting record list.
    export_csv(file_path):
        Exports the pipetting records to a CSV file at the specified file path.
    """
    
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
    """
    Transfers a specified volume of liquid from a source well to a destination well
    using a pipette, with optional mixing steps before and after the transfer, and
    records the pipetting actions with timestamps.
    Args:
        pipette (protocol_api.InstrumentContext): The pipette to perform the transfer.
        volume (float): The volume of liquid to transfer.
        source (protocol_api.labware.Well): The source well to transfer liquid from.
        dest (protocol_api.labware.Well): The destination well to transfer liquid to.
        pipetting_record (PipettingRecordList): The record object to log pipetting actions.
        remark (str, optional): A remark or note to associate with the pipetting action. Defaults to ''.
        mix_before (tuple, optional): A tuple specifying the number of mix cycles and volume 
            to mix in the source well before the transfer. Defaults to None.
        mix_after (tuple, optional): A tuple specifying the number of mix cycles and volume 
            to mix in the destination well after the transfer. Defaults to None.
        drop_tip (bool, optional): Whether to drop the pipette tip after the transfer. Defaults to True.
    Raises:
        None
    Notes:
        - If the pipette does not have a tip, it will automatically pick one up.
        - Mixing steps are optional and can be specified using the `mix_before` and `mix_after` arguments.
        - The pipetting actions are logged in the `pipetting_record` object for tracking purposes.
    """
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


def agar_height_from_weight(source_csv, deck_num):
    """
    Calculate the height of agar in a dish or plate based on its weight.
    This function reads a CSV file containing weight data and calculates the 
    height of agar in either a microplate or a petri dish. The calculation 
    considers the type of plate, its dimensions, and the density of agar.
    Args:
        source_csv (str): Path to the CSV file containing weight data. The CSV 
                          should have no header and contain relevant weight 
                          information.
        deck_num (int): Deck number indicating the position of the data in the 
                        CSV file. Use 6 for the first row and 9 for the second row.
    Returns:
        float: The calculated height of the agar in millimeters (mm). For petri 
               dishes, the dish thickness is added to the calculated height.
    Notes:
        - For microplates:
            - Dish weight is assumed to be 35.28 g.
            - Agar density is 0.00102 g/mm³.
            - Plate bottom area is 9856.785 mm².
        - For petri dishes:
            - Dish weight is assumed to be 17.88 g.
            - Dish thickness is 2 mm.
            - Agar density is 0.00102 g/mm³.
            - Dish diameter is 86 mm.
    """
    
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
