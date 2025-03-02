from opentrons import protocol_api
from commons import (
    transfer_with_timestamp,
    change_pipette_speed,
    PipettingRecordList,
    load_definition_from_json,
    agar_height_from_weight, 
    volume_to_height_50mL_tube,
    list_2d_flatten
)
import pandas as pd
# END_OF_IMPORTS #############################################################

metadata = {
    'apiLevel': '2.13',
    'protocolName': 'yeast RNA-seq (first-half)',
    'description': '''RNA-seq of yeast in iLab''',
    'author': 'Ryosuke Matsuzawa'
}


def run(protocol: protocol_api.ProtocolContext):

    # hyperparameters #########################################################

    USE_SINGLE_TECHREP = True
    BIOREP1_POSITION = 'A4'
    BIOREP2_POSITION = 'B4'
    WATER_POSITION = 'B3'
    TIP300_FIRST_TIP = 'E2'
    TIP20_FIRST_TIP = 'A1'
    # flow rate: [speed1, speed2, speed3, speed4]
    PIPETTE_FLOW_RATE = [50, 130, 210, 290]
    # PIPETTE_FLOW_RATE = [290, 290, 290, 290]
    # MIX_VOLUME = 200
    # MIX_TIMES = 3

    # load dilution recipe
    dilution_recipe = pd.read_csv(
        './liquid_volume.csv',
        sep=',',
        index_col=0
    )
    water_volume = dilution_recipe.iloc[:, 0].values.tolist()
    yeast_volume = dilution_recipe.iloc[:, 1].values.tolist()

    if USE_SINGLE_TECHREP:
        # concat top 8 rows 3 times
        water_volume = water_volume[:8]
        yeast_volume = yeast_volume[:8]
        water_volume = water_volume + water_volume + water_volume
        yeast_volume = yeast_volume + yeast_volume + yeast_volume
    else:
        # use top 24 rows
        water_volume = water_volume[:24]
        yeast_volume = yeast_volume[:24]

    # Config about spotting
    SPOT_SPEED = 50
    SPOT_VOLUME = 2
    # SPOT_HEIGHT = agar_height_from_weight('./agar_plate_weight.csv')
    SPOT_HEIGHT_6 = agar_height_from_weight('./agar_plate_weight.csv', 6)
    SPOT_HEIGHT_9 = agar_height_from_weight('./agar_plate_weight.csv', 9)

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
    p300.starting_tip = tip300.well(TIP300_FIRST_TIP)
    p20.starting_tip = tip20.well(TIP20_FIRST_TIP)

    pipetting_record = PipettingRecordList()

    # Alias to each source
    source_biorep1 = tube_source[BIOREP1_POSITION]
    source_biorep2 = tube_source[BIOREP2_POSITION]
    water = tube_source[WATER_POSITION]

    # Order of Spotting well in microplate
    SPOT_LABEL = [
        'C4', 'D4', 'E4', 'F4',
        'C5', 'D5', 'E5', 'F5',
        'C6', 'D6', 'E6', 'F6',
        'C7', 'D7', 'E7', 'F7',
        'C8', 'D8', 'E8', 'F8',
        'C9', 'D9', 'E9', 'F9'
    ]

    # Order of Edge-Spotting well in microplate
    EDGE_LABEL = [
        'B3', 'C3', 'D3', 'E3', 'F3', 'G3',
        'B4', 'G4', 
        'B5', 'G5',
        'B6', 'G6',
        'B7', 'G7',
        'B8', 'G8',
        'B9', 'G9',
        'B10', 'C10', 'D10', 'E10', 'F10', 'G10'
    ]
    # =========================================================================

    # Order of Spotting
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
    aspirated_water = 0
    # for i in range(24):
    #     surface_height_decrese = volume_to_height_50mL_tube(
    #         water_volume[i]
    #     )
    #     transfer_with_timestamp(
    #         p300,
    #         water_volume[i],
    #         water.bottom(z=95 - aspirated_water),
    #         tube_RNAseq.wells()[i],
    #         pipetting_record,
    #         drop_tip=False,
    #         remark="Add water to RNA-seq tube"
    #     )
    #     aspirated_water += surface_height_decrese

    # Add water to 96 well plate

    if USE_SINGLE_TECHREP:
        end_of_used_well = 8
    else:
        end_of_used_well = 24
    # for i in range(24):
    #     surface_height_decrese = volume_to_height_50mL_tube(
    #         160
    #     )

    #     transfer_with_timestamp(
    #         p300,
    #         200, #180,
    #         water.bottom(z=95 - aspirated_water),
    #         plate_96well.wells()[end_of_used_well + i],
    #         pipetting_record,
    #         drop_tip=False,
    #         remark="Add water to 96 well plate"
    #     )
    #     aspirated_water += surface_height_decrese
    # p300.drop_tip()

    # dilute yeast for rnaseq, for spot assay and spotting
    for index, biorep, speed_index, remark in dilution_order:
        p20.pick_up_tip()
        # change_pipette_speed(
        #     p300,
        #     PIPETTE_FLOW_RATE[speed_index],
        #     PIPETTE_FLOW_RATE[speed_index]
        # )
        # p300.pick_up_tip()
        # p20.pick_up_tip()

        # # Dilute yeast for RNA-seq
        # transfer_with_timestamp(
        #     p300,
        #     yeast_volume[index],
        #     biorep,
        #     # yeast_source_list[yeast_index],
        #     tube_RNAseq.wells()[index].bottom(z=15),
        #     pipetting_record=pipetting_record,
        #     mix_before=(3, 300),
        #     mix_after=(3, 300),
        #     drop_tip=False,
        #     remark=remark + "_RNAseq"
        # )

        # # Dilute yeast for spot assay
        # transfer_with_timestamp(
        #     p300,
        #     20,
        #     tube_RNAseq.wells()[index].bottom(z=15),
        #     plate_96well.wells()[end_of_used_well + index],
        #     pipetting_record=pipetting_record,
        #     mix_after=(3, 120),
        #     drop_tip=False,
        #     remark=remark + "_96well"
        # )

        # # Spotting on Deck6 agar plate
        # change_pipette_speed(p20, SPOT_SPEED, SPOT_SPEED)
        # agar_plate_row_flatten = list_2d_flatten(agar_plate_6.rows())
        # transfer_with_timestamp(
        #     p20,
        #     SPOT_VOLUME,
        #     plate_96well.wells()[end_of_used_well + index],
        #     agar_plate_row_flatten[index].bottom(z=SPOT_HEIGHT_6),
        #     pipetting_record=pipetting_record,
        #     drop_tip=False,
        #     remark=remark + "_spotting"
        # )
        
        
        # Spotting on Deck9 agar plate
        transfer_with_timestamp(
            p20,
            SPOT_VOLUME,
            plate_96well.wells()[end_of_used_well + index],
            agar_plate_9.wells_by_name()[SPOT_LABEL[index]].bottom(z=SPOT_HEIGHT_9),
            pipetting_record=pipetting_record,
            drop_tip=True,
            remark=remark + "_spotting"
        )

        # p300.drop_tip()
        # p20.drop_tip()

    for index, biorep, speed_index, remark in dilution_order:
        p20.pick_up_tip()
        
        # Spotting on Deck9 agar plate
        transfer_with_timestamp(
            p20,
            SPOT_VOLUME,
            plate_96well.wells()[end_of_used_well + index],
            agar_plate_9.wells_by_name()[EDGE_LABEL[index]].bottom(z=SPOT_HEIGHT_9),
            pipetting_record=pipetting_record,
            drop_tip=True,
            remark=remark + "_edge-spotting"
        )
        
            
    # Output record
    pipetting_record.export_csv('pipetting_record_second.csv')
