from opentrons import protocol_api

metadata = {
    'apiLevel': '2.13',
    'protocolName': 'homing',
    'description': '''homing before any runs''',
    'author': 'Shodai Taguchi'
}

def run(protocol: protocol_api.ProtocolContext):
    protocol.home()