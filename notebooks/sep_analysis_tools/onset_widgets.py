"""
A library to run the interactive user interface in SEP event onset determination notebooks.

@ Author: Christian Palmroos <chospa@utu.fi>
@Last updated: 2022-09-30
"""


from importlib.resources import path
import ipywidgets as widgets

# a list of available spacecraft:
# list_of_sc = ["STEREO-A", "STEREO-B", "Solar Orbiter", "Bepicolombo", "SOHO"]
list_of_sc = ["STEREO-A", "STEREO-B", "Solar Orbiter", "SOHO"]

stereo_instr = ["SEPT", "HET"] #["LET", "SEPT", "HET"]
solo_instr = ["EPT", "HET"]
bepi_instr = ["SIXS-P"]
soho_instr = ["ERNE-HED", "EPHIN"]

sensor_dict = {
    "STEREO-A" : stereo_instr,
    "STEREO-B" : stereo_instr,
    "Solar Orbiter" : solo_instr,
    "Bepicolombo" : bepi_instr,
    "SOHO" : soho_instr
}

view_dict = {
    ("STEREO-A", "SEPT") : ["sun", "asun", "north", "south"],
    ("STEREO-B", "SEPT") : ["sun", "asun", "north", "south"],
    ("Solar Orbiter", "EPT") : ["sun", "asun", "north", "south"],
    ("Solar Orbiter", "HET") : ["sun", "asun", "north", "south"],
    ("Bepicolombo", "SIXS-P") : [0, 1, 2, 3, 4]
}

species_dict = {
    ("STEREO-A", "LET") : ['protons', 'electrons'],
    ("STEREO-A", "SEPT") : ['ions', 'electrons'],
    ("STEREO-A", "HET") : ['protons', 'electrons'],
    ("STEREO-B", "LET") : ['protons', 'electrons'],
    ("STEREO-B", "SEPT") : ['ions', 'electrons'],
    ("STEREO-B", "HET") : ['protons', 'electrons'],
    ("Solar Orbiter", "EPT") : ['ions', 'electrons'],
    ("Solar Orbiter", "HET") : ['protons', 'electrons'],
    ("Bepicolombo", "SIXS-P") : ['protons', 'electrons'],
    ("SOHO", "ERNE-HED") : ['protons'],
    ("SOHO", "EPHIN") : ['electrons']
}

radio_dict = {
    "None" : None,
    "STEREO-A" : ("ahead", "STEREO-A"),
    "STEREO-B" : ("behind", "STEREO-B"),
    "WIND (Coming soon!)" : ("wind", "WIND")
}

# Drop-downs for dynamic particle spectrum:
spacecraft_drop = widgets.Dropdown(
                                options = list_of_sc,
                                description = "Spacecraft:",
                                disabled = False,
                                )

sensor_drop = widgets.Dropdown(
                                options = sensor_dict[spacecraft_drop.value],
                                description = "Sensor:",
                                disabled = False,
                                )

view_drop = widgets.Dropdown(
                                options = view_dict[(spacecraft_drop.value, sensor_drop.value)],
                                description = "Viewing:",
                                disabled = False
                                )

species_drop = widgets.Dropdown(
                                options = species_dict[(spacecraft_drop.value, sensor_drop.value)],
                                description = "Species:",
                                disabled = False,
                                )


# A button to enable radio spectrum (Leave this out for now, sincde it doesn't work in the server as of 2022-09-30)
radio_button = widgets.Checkbox(
                                value=False,
                                description='Radio Spectrum',
                                disabled=True,
                                indent=False
                                )

# The drop-drown for radio options
radio_drop_style = {'description_width': 'initial'}
radio_drop = widgets.Dropdown(
                                options = radio_dict,
                                value = None,
                                description = "Plot radio spectrum for:",
                                disabled = False,
                                style=radio_drop_style
                              )

# this function updates the options in sensor_drop menu
def update_sensor_options(val):
    sensor_drop.options = sensor_dict[spacecraft_drop.value]


# updates the options and availability of view_drop menu
def update_view_options(val):
    try:
        view_drop.disabled = False
        view_drop.options = view_dict[(spacecraft_drop.value, sensor_drop.value)]
        view_drop.value = view_drop.options[0]
    except KeyError:
        view_drop.disabled = True
        view_drop.value = None


def update_species_options(val):
    try:
        species_drop.options = species_dict[(spacecraft_drop.value, sensor_drop.value)]
    except KeyError:
        pass


def update_radio_options(val):
    radio_drop.disabled = not radio_button.value
    if radio_drop.disabled:
        radio_drop.value = None
    else:
        radio_drop.value = radio_drop.options[0]


def confirm_input(event_date : int, data_path : str, plot_path : str):

    print("You've chosen the following options:")
    print(f"Spacecraft: {spacecraft_drop.value}")
    print(f"Sensor: {sensor_drop.value}")
    print(f"Species: {species_drop.value}")
    print(f"Viewing: {view_drop.value}")
    print(f"Event_date: {event_date}")
    print(f"Data_path: {data_path}")
    print(f"Plot_path: {plot_path}")

    if spacecraft_drop.value == "Solar Orbiter":
        spacecraft_drop_value = "solo"
    elif spacecraft_drop.value == "STEREO-A":
        spacecraft_drop_value = "sta"
    elif spacecraft_drop.value == "STEREO-B":
        spacecraft_drop_value = "stb"
    else:
        spacecraft_drop_value = spacecraft_drop.value
    
    if sensor_drop.value in ["ERNE-HED"]:
        sensor_drop_value = "ERNE"
    else:
        sensor_drop_value = sensor_drop.value
    
    if species_drop.value == "protons":
        species_drop_value = 'p'
    else:
        species_drop_value = 'e'

    # this is to be fed into Event class as input
    global input_dict

    input_dict = {
        "Spacecraft" : spacecraft_drop_value,
        "Sensor" : sensor_drop_value,
        "Species" : species_drop_value,
        "Viewing" : view_drop.value,
        "Event_date" : event_date,
        "Data_path" : data_path,
        "Plot_path" : plot_path
    }

# makes spacecraft_drop run these functions every time it is accessed by user
spacecraft_drop.observe(update_sensor_options)
spacecraft_drop.observe(update_view_options)
sensor_drop.observe(update_view_options)

# does the same but for sensor menu
spacecraft_drop.observe(update_species_options)
sensor_drop.observe(update_species_options)

# also observe the radio menu
# radio_button.observe(update_radio_options)
