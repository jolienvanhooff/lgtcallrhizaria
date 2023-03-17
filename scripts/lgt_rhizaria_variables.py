#!/usr/bin/env python3

# Module used to import various variables - such as dictionaries - often used in LgtRizaria project

eventcolors = {
    "lateral": '#DC143C', # Crimson
    "vertical": '#00FFFF', # Aqua
    "duplication": '#800000', # Maroon
    "losses": '#4B0082', # Indigo
}

taxa_extant = ['ParGloRCC365', 'LotGloLEX01', 'Lot_CCMP622', 'LotAmoCCMP2058', 'Gym_CCMP2014', 'NorSphBC52', 'BigNatCCMP2755', 'BigLonCCMP242', 'ChlRep', 'GutVulBF0784', 'AurSolAF21', 'PauChr', 'MatSpe', 'AboPro', 'BreSpeRG2016a', 'PlaBraE3', 'SpoSub', 'LepVor', 'GroSph', 'MikMac', 'LapGus', 'RetFil', 'AmmSpe', 'SorSpe', 'ElpMar', 'RosSpe', 'GloSpeGF15', 'StiZan', 'LitSet']

non_terminal_and_genomic_taxa = ['Rhizaria', 'Cercozoa', 'Retaria', 'Cercozoa1', 'Endomyxa', 'Foraminifera', 'Radiolaria', 'Reticulofilosa', 'Monadofilosa', 'Endomyxa2', 'Globothalamea', 'Chlorarachniophyceae2', 'Chlorarachniophyceae1', 'Monadofilosa1', 'Endomyxa3', 'Ascetosporea', 'Rotaliida', 'Chlorarachniophyceae3', 'Chlorarachniophyceae4', 'Lotharella', 'Monadofilosa4', 'Monadofilosa2', 'Plasmodiophoridae', 'Rotaliida1', 'Bigelowiella', 'Chlorarachniophyceae5', 'Monadofilosa3', 'Rotalioidea', 'PlaBraE3', 'BigNatCCMP2755']

taxa_ancestral = ['Rhizaria', 'Cercozoa', 'Retaria', 'Cercozoa1', 'Endomyxa', 'Foraminifera', 'Radiolaria', 'Reticulofilosa', 'Monadofilosa', 'Endomyxa2', 'Globothalamea', 'Chlorarachniophyceae2', 'Chlorarachniophyceae1', 'Monadofilosa1', 'Endomyxa3', 'Ascetosporea', 'Rotaliida', 'Chlorarachniophyceae3', 'Chlorarachniophyceae4', 'Lotharella', 'Monadofilosa4', 'Monadofilosa2', 'Plasmodiophoridae', 'Rotaliida1', 'Bigelowiella', 'Chlorarachniophyceae5', 'Monadofilosa3', 'Rotalioidea']

taxa_extant_ordered = ['BigNatCCMP2755', 'BigLonCCMP242', 'NorSphBC52', 'ChlRep', 'Gym_CCMP2014', 'LotAmoCCMP2058', 'ParGloRCC365', 'LotGloLEX01', 'Lot_CCMP622', 'AurSolAF21', 'BreSpeRG2016a', 'PauChr', 'MatSpe', 'AboPro', 'GutVulBF0784', 'PlaBraE3', 'SpoSub', 'LepVor', 'GroSph', 'MikMac', 'LapGus', 'AmmSpe', 'ElpMar', 'RosSpe', 'GloSpeGF15', 'SorSpe', 'RetFil', 'StiZan', 'LitSet']

taxa_all_ordered = ["Rhizaria", "Cercozoa", "Cercozoa1", "Reticulofilosa", "Chlorarachniophyceae2", "Chlorarachniophyceae3", "Bigelowiella", "BigNatCCMP2755", "BigLonCCMP242", "NorSphBC52", "ChlRep", "Chlorarachniophyceae1", "Chlorarachniophyceae4", "Chlorarachniophyceae5", "Gym_CCMP2014", "LotAmoCCMP2058", "ParGloRCC365", "Lotharella", "LotGloLEX01", "Lot_CCMP622", "Monadofilosa", "Monadofilosa1", "Monadofilosa4", "AurSolAF21", "BreSpeRG2016a", "Monadofilosa2", "PauChr", "Monadofilosa3", "MatSpe", "AboPro", "GutVulBF0784", "Endomyxa", "Endomyxa2", "Endomyxa3", "Plasmodiophoridae", "PlaBraE3", "SpoSub", "LepVor", "Ascetosporea", "GroSph", "MikMac", "LapGus", "Retaria", "Foraminifera", "Globothalamea", "Rotaliida", "Rotaliida1", "Rotalioidea", "AmmSpe", "ElpMar", "RosSpe", "GloSpeGF15", "SorSpe", "RetFil", "Radiolaria", "StiZan", "LitSet"]

strainids = {
    "ParGloRCC365" : "Partenskyella glossopodia Strain RCC365",
    "LotGloLEX01" : "Lotharella globosa LEX01",
    "Lot_CCMP622" : "Lotharella oceanica CCMP622",
    "LotAmoCCMP2058" : "Amorphochlora amoebiformis Strain CCMP2058",
    "Gym_CCMP2014" : "Gymnochlora CCMP2014",
    "NorSphBC52" : "Norrisiella sphaerica Strain BC52",
    "BigNatCCMP2755" : "Bigelowiella natans CCMP2755",
    "BigLonCCMP242" : "Bigelowiella longifila Strain CCMP242",
    "ChlRep" : "Chlorarachnion reptans",
    "GutVulBF0784" : "Guttulinopsis vulgaris BF07-8-4",
    "AurSolAF21" : "Aurigamonas solis AF-21",
    "PauChr" : "Paulinella chromatophora",
    "MatSpe" : "Mataza sp. / unidentified eukaryote D1",
    "AboPro" : "Abollifer prolabens",
    "BreSpeRG2016a" : "Brevimastigomonas sp. RG-2016a",
    "PlaBraE3" : "Plasmodiophora brassicae isolate e3",
    "SpoSub" : "Spongospora subterranea",
    "LepVor" : "Leptophrys vorax",
    "GroSph" : "Gromia sphaerica",
    "MikMac" : "Mikrocytos mackini",
    "LapGus" : "Lapot gusevi",
    "RetFil" : "Reticulomyxa filosa",
    "AmmSpe" : "Ammonia sp.",
    "SorSpe" : "Sorites sp.",
    "ElpMar" : "Elphidium margaritaceum",
    "RosSpe" : "Rosalina sp.",
    "GloSpeGF15" : "Globobulimina sp. GF15",
    "StiZan" : "Sticholonche zanclea",
    "LitSet" : "Lithomelissa setosa"
}

speciesids = {
    "ParGloRCC365" : "Partenskyella glossopodia",
    "LotGloLEX01" : "Lotharella globosa",
    "Lot_CCMP622" : "Lotharella oceanica",
    "LotAmoCCMP2058" : "Amorphochlora amoebiformis",
    "Gym_CCMP2014" : "Gymnochlora sp.",
    "NorSphBC52" : "Norrisiella sphaerica",
    "BigNatCCMP2755" : "Bigelowiella natans",
    "BigLonCCMP242" : "Bigelowiella longifila",
    "ChlRep" : "Chlorarachnion reptans",
    "GutVulBF0784" : "Guttulinopsis vulgaris",
    "AurSolAF21" : "Aurigamonas solis",
    "PauChr" : "Paulinella chromatophora",
    "MatSpe" : "Mataza sp.",
    "AboPro" : "Abollifer prolabens",
    "BreSpeRG2016a" : "Brevimastigomonas sp.",
    "PlaBraE3" : "Plasmodiophora brassicae",
    "SpoSub" : "Spongospora subterranea",
    "LepVor" : "Leptophrys vorax",
    "GroSph" : "Gromia sphaerica",
    "MikMac" : "Mikrocytos mackini",
    "LapGus" : "Lapot gusevi",
    "RetFil" : "Reticulomyxa filosa",
    "AmmSpe" : "Ammonia sp.",
    "SorSpe" : "Sorites sp.",
    "ElpMar" : "Elphidium margaritaceum",
    "RosSpe" : "Rosalina sp.",
    "GloSpeGF15" : "Globobulimina sp.",
    "StiZan" : "Sticholonche zanclea",
    "LitSet" : "Lithomelissa setosa"
}
