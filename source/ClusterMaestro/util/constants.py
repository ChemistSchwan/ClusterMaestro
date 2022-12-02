"""Provides various constants for internal handling og e.g. stomic distances...
"""
colorList_1 = [
    "#005b9a",
    "#15b434",
    "#c33232",
    "#972b9a",
    "#46c1c0",
    "#c07f2b",
    "#005b9a",
    "#15b434",
]

colorList = colorList_1

colorList_2 = [
    "#005b9a",
    "#15b434",
    "#c33232",
    "#972b9a",
    "#46c1c0",
    "#c07f2b",
    "#007ea0",
    "#941651",
    "#ffb000",
    "#696966",
    "#60ff05",
    "#ff5c1c",
    "#c802f9",
    "#0f6c75",
    "#fcf80c",
    "#fc0c78",
    "#7c7c7c",
]

colorList_3 = [
    "#7e1000",
    "#ff0000",
    "#ff8585",
    "#001b7e",
    "#002cff",
    "#85a0ff",
    "#007e03",
    "#05cd00",
    "#49f833",
    "#7e6800",
    "#ffc800",
    "#e7dc5f",
]

ptable = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Uut",
    114: "Fl",
    115: "Uup",
    116: "Lv",
    117: "Uus",
    118: "Uuo",
}

atomMasses = {
    "H": 1.0079,
    "He": 4.0026,
    "Li": 6.941,
    "Be": 9.0122,
    "B": 10.811,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9984,
    "F": 18.9984,
    "Ne": 20.1797,
    "Na": 22.990,
    "Mg": 24.305,
    "Al": 26.981,
    "Si": 28.085,
    "P": 30.974,
    "S": 32.067,
    "Cl": 35.451,
    "Ar": 39.948,
    "K": 39.098,
    "Ca": 40.078,
    "Sc": 44.956,
    "Ti": 47.867,
    "V": 50.941,
    "Cr": 51.996,
    "Mn": 54.938,
    "Fe": 55.945,
    "Co": 58.933,
    "Ni": 58.693,
    "Cu": 63.546,
    "Zn": 65.380,
    "Ga": 69.723,
    "Ge": 72.631,
    "As": 74.922,
    "Se": 78.972,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.468,
    "Sr": 87.620,
    "Y": 88.906,
    "Zr": 91.224,
    "Nb": 92.906,
    "Mo": 95.950,
    "Tc": 97.000,
    "Ru": 101.07,
    "Rh": 102.91,
    "Pd": 106.42,
    "Ag": 107.868,
    "Cd": 112.414,
    "In": 114.818,
    "Sn": 118.711,
    "Sb": 121.760,
    "Te": 127.760,
    "I": 126.904,
    "Xe": 131.294,
    "Cs": 132.905,
    "Ba": 137.328,
    "La": 138.905,
    "Hf": 178.486,
    "Ta": 180.948,
    "W": 183.841,
    "Re": 186.207,
    "Os": 170.233,
    "Ir": 195.217,
    "PT": 195.084,
    "Au": 196.967,
    "Hg": 200.592,
    "Tl": 204.38,
    "Pb": 207.21,
    "Bi": 208.980,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Ce": 140.116,
    "Pr": 140.907,
    "Nd": 144.242,
    "Pm": 145,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.925,
    "Dy": 162.500,
    "Ho": 164.930,
    "Er": 167.259,
    "Tm": 168.934,
    "Yb": 173.045,
    "Lu": 174.967,
}

#               '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,    '<++>': <++>,  '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,
#               '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,    '<++>': <++>,  '<++>': <++>,   '<++>': <++>,   '<++>': <++>,   '<++>': <++>,

inv_ptable = {v: k for k, v in ptable.items()}

vdw_radii = {1: 1.0, 6: 1.8, 14: 2.25, 16: 2.6, 32: 2.35, 50: 2.6}


def get_scaling(new_atom_spec):
    rings, bridges, corners = {
        "Sn": (1.27, 1.587, 1.491),
        "Ge": (1.24, 1.5, 1.44),
        "Si": (1.23, 1.42, 1.38),
    }.get(new_atom_spec, (1.0, 1.0, 1.0))
    return rings, bridges, corners


rad_dict = {
    "H": 0.2,
    "He": 0.4,
    "Li": 0.4,
    "Be": 0.4,
    "B": 0.4,
    "C": 0.32,
    "N": 0.4,
    "O": 0.4,
    "F": 0.4,
    "Ne": 0.4,
    "Na": 0.4,
    "Mg": 0.4,
    "Al": 0.4,
    "Si": 0.4,
    "P": 0.4,
    "S": 0.4,
    "Cl": 0.4,
    "Ar": 0.4,
    "K ": 0.4,
    "Ca": 0.4,
    "Sc": 0.4,
    "Ti": 0.4,
    "V": 0.4,
    "Cr": 0.4,
    "Mn": 0.4,
    "Fe": 0.4,
    "Co": 0.4,
    "Ni": 0.4,
    "Cu": 0.4,
    "Zn": 0.4,
    "Ga": 0.4,
    "Ge": 0.4,
    "As": 0.4,
    "Se": 0.4,
    "Br": 0.4,
    "Kr": 0.4,
    "Rb": 0.4,
    "Sr": 0.4,
    "Y": 0.4,
    "Zr": 0.4,
    "Nb": 0.4,
    "Mo": 0.4,
    "Tc": 0.4,
    "Ru": 0.4,
    "Rh": 0.4,
    "Pd": 0.4,
    "Ag": 0.4,
    "Cd": 0.4,
    "In": 0.4,
    "Sn": 0.4,
    "Sb": 0.4,
    "Te": 0.4,
    "I": 0.4,
    "Xe": 0.4,
    "Cs": 0.4,
    "Ba": 0.4,
    "La": 0.4,
    "Ce": 0.4,
    "Pr": 0.4,
    "Nd": 0.4,
    "Pm": 0.4,
    "Sm": 0.4,
    "Eu": 0.4,
    "Gd": 0.4,
    "Tb": 0.4,
    "Dy": 0.4,
    "Ho": 0.4,
    "Er": 0.4,
    "Tm": 0.4,
    "Yb": 0.4,
    "Lu": 0.4,
    "Hf": 0.4,
    "Ta": 0.4,
    "W": 0.4,
    "Re": 0.4,
    "Os": 0.4,
    "Ir": 0.4,
    "Pt": 0.4,
    "Au": 0.4,
    "Hg": 0.4,
    "Tl": 0.4,
    "Pb": 0.4,
    "Bi": 0.4,
    "Po": 0.4,
    "At": 0.4,
    "Rn": 0.4,
    "Fr": 0.4,
    "Ra": 0.4,
    "Ac": 0.4,
    "Th": 0.4,
    "Pa": 0.4,
    "U": 0.4,
    "Np": 0.4,
    "Pu": 0.4,
    "Am": 0.4,
    "Cm": 0.4,
    "Bk": 0.4,
    "Cf": 0.4,
    "Es": 0.4,
    "Fm": 0.4,
    "Md": 0.4,
    "No": 0.4,
    "Lr": 0.4,
    "Rf": 0.4,
    "Db": 0.4,
    "Sg": 0.4,
    "Bh": 0.4,
    "Hs": 0.4,
    "Mt": 0.4,
}


rgb_dict = {
    "H": [217, 255, 255],
    "He": [255, 255, 160],
    "Li": [204, 128, 255],
    "Be": [194, 255, 0],
    "B": [255, 181, 181],
    "C": [144, 144, 144],
    "N": [48, 80, 248],
    "O": [255, 13, 13],
    "F": [144, 224, 80],
    "Ne": [179, 227, 245],
    "Na": [171, 92, 242],
    "Mg": [138, 255, 0],
    "Al": [191, 166, 166],
    "Si": [240, 200, 160],
    "P": [255, 128, 0],
    "S": [255, 255, 48],
    "Cl": [31, 240, 31],
    "Ar": [128, 209, 227],
    "K ": [143, 64, 212],
    "Ca": [61, 255, 0],
    "Sc": [230, 230, 230],
    "Ti": [191, 194, 199],
    "V": [166, 166, 171],
    "Cr": [138, 153, 199],
    "Mn": [156, 122, 199],
    "Fe": [224, 102, 51],
    "Co": [240, 144, 160],
    "Ni": [80, 208, 80],
    "Cu": [200, 128, 51],
    "Zn": [125, 128, 176],
    "Ga": [194, 143, 143],
    "Ge": [102, 143, 143],
    "As": [189, 128, 227],
    "Se": [255, 161, 0],
    "Br": [166, 41, 41],
    "Kr": [92, 184, 209],
    "Rb": [112, 46, 176],
    "Sr": [0, 255, 0],
    "Y": [148, 255, 255],
    "Zr": [148, 224, 224],
    "Nb": [115, 194, 201],
    "Mo": [84, 181, 181],
    "Tc": [59, 158, 158],
    "Ru": [36, 143, 143],
    "Rh": [10, 125, 140],
    "Pd": [0, 105, 133],
    "Ag": [192, 192, 192],
    "Cd": [255, 217, 143],
    "In": [166, 117, 115],
    "Sn": [102, 128, 128],
    "Sb": [158, 99, 181],
    "Te": [212, 122, 0],
    "I": [148, 0, 148],
    "Xe": [66, 158, 176],
    "Cs": [87, 23, 143],
    "Ba": [0, 201, 0],
    "La": [112, 212, 255],
    "Ce": [255, 255, 199],
    "Pr": [217, 255, 199],
    "Nd": [199, 255, 199],
    "Pm": [163, 255, 199],
    "Sm": [143, 255, 199],
    "Eu": [97, 255, 199],
    "Gd": [69, 255, 199],
    "Tb": [48, 255, 199],
    "Dy": [31, 255, 199],
    "Ho": [0, 255, 156],
    "Er": [0, 230, 117],
    "Tm": [0, 212, 82],
    "Yb": [0, 191, 56],
    "Lu": [0, 171, 36],
    "Hf": [77, 194, 255],
    "Ta": [77, 166, 255],
    "W": [33, 148, 214],
    "Re": [38, 125, 171],
    "Os": [38, 102, 150],
    "Ir": [23, 84, 135],
    "Pt": [208, 208, 224],
    "Au": [255, 209, 35],
    "Hg": [184, 184, 208],
    "Tl": [166, 84, 77],
    "Pb": [87, 89, 97],
    "Bi": [158, 79, 181],
    "Po": [171, 92, 0],
    "At": [117, 79, 69],
    "Rn": [66, 130, 150],
    "Fr": [66, 0, 102],
    "Ra": [0, 125, 0],
    "Ac": [112, 171, 250],
    "Th": [0, 186, 255],
    "Pa": [0, 161, 255],
    "U": [0, 143, 255],
    "Np": [0, 128, 255],
    "Pu": [0, 107, 255],
    "Am": [84, 92, 242],
    "Cm": [120, 92, 227],
    "Bk": [138, 79, 227],
    "Cf": [161, 54, 212],
    "Es": [179, 31, 212],
    "Fm": [179, 31, 186],
    "Md": [179, 13, 166],
    "No": [189, 13, 135],
    "Lr": [199, 0, 102],
    "Rf": [204, 0, 89],
    "Db": [209, 0, 79],
    "Sg": [217, 0, 69],
    "Bh": [224, 0, 56],
    "Hs": [230, 0, 46],
    "Mt": [235, 0, 38],
}

bond_dist_dict = {
    "H": 1.0,  # default
    "He": 2.0,  # default
    "Li": 2.0,  # default
    "Be": 2.0,  # default
    "B": 2.0,  # default
    "C": 1.5,
    "N": 2.0,  # default
    "O": 2.0,  # default
    "F": 2.0,  # default
    "Ne": 2.0,  # default
    "Na": 2.0,  # default
    "Mg": 2.0,  # default
    "Al": 2.0,  # default
    "Si": 2.25,
    "P": 2.0,  # default
    "S": 2.1,
    "Cl": 2.0,  # default
    "Ar": 2.0,  # default
    "K ": 2.0,  # default
    "Ca": 2.0,  # default
    "Sc": 2.0,  # default
    "Ti": 2.0,  # default
    "V": 2.0,  # default
    "Cr": 2.0,  # default
    "Mn": 2.0,  # default
    "Fe": 2.0,  # default
    "Co": 2.0,  # default
    "Ni": 2.0,  # default
    "Cu": 2.0,  # default
    "Zn": 2.0,  # default
    "Ga": 2.0,  # default
    "Ge": 2.35,
    "As": 2.0,  # default
    "Se": 2.45,  # default
    "Br": 2.0,  # default
    "Kr": 2.0,  # default
    "Rb": 2.0,  # default
    "Sr": 2.0,  # default
    "Y": 2.0,  # default
    "Zr": 2.0,  # default
    "Nb": 2.0,  # default
    "Mo": 2.0,  # default
    "Tc": 2.0,  # default
    "Ru": 2.0,  # default
    "Rh": 2.0,  # default
    "Pd": 2.0,  # default
    "Ag": 2.0,  # default
    "Cd": 2.0,  # default
    "In": 2.0,  # default
    "Sn": 2.7,  # Custom elongated for connectivity matrix
    "Sb": 2.0,  # default
    "Te": 2.9,
    "I": 2.0,  # default
    "Xe": 2.0,  # default
    "Cs": 2.0,  # default
    "Ba": 2.0,  # default
    "La": 2.0,  # default
    "Ce": 2.0,  # default
    "Pr": 2.0,  # default
    "Nd": 2.0,  # default
    "Pm": 2.0,  # default
    "Sm": 2.0,  # default
    "Eu": 2.0,  # default
    "Gd": 2.0,  # default
    "Tb": 2.0,  # default
    "Dy": 2.0,  # default
    "Ho": 2.0,  # default
    "Er": 2.0,  # default
    "Tm": 2.0,  # default
    "Yb": 2.0,  # default
    "Lu": 2.0,  # default
    "Hf": 2.0,  # default
    "Ta": 2.0,  # default
    "W": 2.0,  # default
    "Re": 2.0,  # default
    "Os": 2.0,  # default
    "Ir": 2.0,  # default
    "Pt": 2.0,  # default
    "Au": 2.0,  # default
    "Hg": 2.0,  # default
    "Tl": 2.0,  # default
    "Pb": 2.0,  # default
    "Bi": 2.0,  # default
    "Po": 2.0,  # default
    "At": 2.0,  # default
    "Rn": 2.0,  # default
    "Fr": 2.0,  # default
    "Ra": 2.0,  # default
    "Ac": 2.0,  # default
    "Th": 2.0,  # default
    "Pa": 2.0,  # default
    "U": 2.0,  # default
    "Np": 2.0,  # default
    "Pu": 2.0,  # default
    "Am": 2.0,  # default
    "Cm": 2.0,  # default
    "Bk": 2.0,  # default
    "Cf": 2.0,  # default
    "Es": 2.0,  # default
    "Fm": 2.0,  # default
    "Md": 2.0,  # default
    "No": 2.0,  # default
    "Lr": 2.0,  # default
    "Rf": 2.0,  # default
    "Db": 2.0,  # default
    "Sg": 2.0,  # default
    "Bh": 2.0,  # default
    "Hs": 2.0,  # default
    "Mt": 2.0,
}  # default

plot_dict = {
    "H": [10, "#0099ff"],  # default
    "He": [10, "#ffffff"],  # default
    "Li": [10, "#ffffff"],  # default
    "Be": [10, "#ffffff"],  # default
    "B": [10, "#ffffff"],  # default
    "C": [30, "#663300"],
    "N": [10, "#ffffff"],  # default
    "O": [35, "#ff0000"],
    "F": [10, "#ffffff"],  # default
    "Ne": [10, "#ffffff"],  # default
    "Na": [10, "#ffffff"],  # default
    "Mg": [10, "#ffffff"],  # default
    "Al": [10, "#ffffff"],  # default
    "Si": [10, "#ffffff"],  # default
    "P": [10, "#ffffff"],  # default
    "S": [10, "#ffffff"],  # default
    "Cl": [10, "#ffffff"],  # default
    "Ar": [10, "#ffffff"],  # default
    "K ": [10, "#ffffff"],  # default
    "Ca": [10, "#ffffff"],  # default
    "Sc": [10, "#ffffff"],  # default
    "Ti": [10, "#ffffff"],  # default
    "V": [10, "#ffffff"],  # default
    "Cr": [10, "#ffffff"],  # default
    "Mn": [10, "#ffffff"],  # default
    "Fe": [10, "#ffffff"],  # default
    "Co": [10, "#ffffff"],  # default
    "Ni": [10, "#ffffff"],  # default
    "Cu": [10, "#ffffff"],  # default
    "Zn": [10, "#ffffff"],  # default
    "Ga": [10, "#ffffff"],  # default
    "Ge": [10, "#ffffff"],  # default
    "As": [10, "#ffffff"],  # default
    "Se": [10, "#ffffff"],  # default
    "Br": [10, "#ffffff"],  # default
    "Kr": [10, "#ffffff"],  # default
    "Rb": [10, "#ffffff"],  # default
    "Sr": [10, "#ffffff"],  # default
    "Y": [10, "#ffffff"],  # default
    "Zr": [10, "#ffffff"],  # default
    "Nb": [10, "#ffffff"],  # default
    "Mo": [10, "#ffffff"],  # default
    "Tc": [10, "#ffffff"],  # default
    "Ru": [10, "#ffffff"],  # default
    "Rh": [10, "#ffffff"],  # default
    "Pd": [10, "#ffffff"],  # default
    "Ag": [10, "#ffffff"],  # default
    "Cd": [10, "#ffffff"],  # default
    "In": [10, "#ffffff"],  # default
    "Sn": [10, "#ffffff"],  # default
    "Sb": [10, "#ffffff"],  # default
    "Te": [10, "#ffffff"],  # default
    "I": [10, "#ffffff"],  # default
    "Xe": [10, "#ffffff"],  # default
    "Cs": [10, "#ffffff"],  # default
    "Ba": [10, "#ffffff"],  # default
    "La": [10, "#ffffff"],  # default
    "Ce": [10, "#ffffff"],  # default
    "Pr": [10, "#ffffff"],  # default
    "Nd": [10, "#ffffff"],  # default
    "Pm": [10, "#ffffff"],  # default
    "Sm": [10, "#ffffff"],  # default
    "Eu": [10, "#ffffff"],  # default
    "Gd": [10, "#ffffff"],  # default
    "Tb": [10, "#ffffff"],  # default
    "Dy": [10, "#ffffff"],  # default
    "Ho": [10, "#ffffff"],  # default
    "Er": [10, "#ffffff"],  # default
    "Tm": [10, "#ffffff"],  # default
    "Yb": [10, "#ffffff"],  # default
    "Lu": [10, "#ffffff"],  # default
    "Hf": [10, "#ffffff"],  # default
    "Ta": [10, "#ffffff"],  # default
    "W": [10, "#ffffff"],  # default
    "Re": [10, "#ffffff"],  # default
    "Os": [10, "#ffffff"],  # default
    "Ir": [10, "#ffffff"],  # default
    "Pt": [10, "#ffffff"],  # default
    "Au": [10, "#ffffff"],  # default
    "Hg": [10, "#ffffff"],  # default
    "Tl": [10, "#ffffff"],  # default
    "Pb": [10, "#ffffff"],  # default
    "Bi": [10, "#ffffff"],  # default
    "Po": [10, "#ffffff"],  # default
    "At": [10, "#ffffff"],  # default
    "Rn": [10, "#ffffff"],  # default
    "Fr": [10, "#ffffff"],  # default
    "Ra": [10, "#ffffff"],  # default
    "Ac": [10, "#ffffff"],  # default
    "Th": [10, "#ffffff"],  # default
    "Pa": [10, "#ffffff"],  # default
    "U": [10, "#ffffff"],  # default
    "Np": [10, "#ffffff"],  # default
    "Pu": [10, "#ffffff"],  # default
    "Am": [10, "#ffffff"],  # default
    "Cm": [10, "#ffffff"],  # default
    "Bk": [10, "#ffffff"],  # default
    "Cf": [10, "#ffffff"],  # default
    "Es": [10, "#ffffff"],  # default
    "Fm": [10, "#ffffff"],  # default
    "Md": [10, "#ffffff"],  # default
    "No": [10, "#ffffff"],  # default
    "Lr": [10, "#ffffff"],  # default
    "Rf": [10, "#ffffff"],  # default
    "Db": [10, "#ffffff"],  # default
    "Sg": [10, "#ffffff"],  # default
    "Bh": [10, "#ffffff"],  # default
    "Hs": [10, "#ffffff"],  # default
    "Mt": [10, "#ffffff"],
}  # default
