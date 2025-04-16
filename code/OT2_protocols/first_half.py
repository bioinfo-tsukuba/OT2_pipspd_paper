from opentrons import protocol_api
from commons import (
    transfer_with_timestamp,
    change_pipette_speed,
    PipettingRecordList,
    load_definition_from_json
)
# END_OF_IMPORTS #############################################################


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
    PIPETTE_FLOW_RATE = [50, 130, 210, 290]
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
