"""
Builds data/metabolite_library.csv — curated annotation library for central carbon / nitrogen metabolism.
IDs were curated by hand; fields left blank when not certain. The app's "Verify / fill IDs online" step checks
and completes them against KEGG REST and PubChem PUG REST from the user's computer.
Pathway codes:
 GLY Glycolysis · PPP Pentose phosphate pathway · TCA TCA cycle · GLN Glutamine metabolism ·
 GLU Glutamate metabolism · AA Amino acid metabolism · FA Fatty acid metabolism · PUR Purine metabolism ·
 PYR Pyrimidine metabolism · OC One-carbon metabolism · RED Redox metabolism · NUC Nucleotide biosynthesis ·
 HEX Hexosamine & nucleotide sugars · PL Phospholipid head groups · UREA Urea cycle & polyamines ·
 COF Cofactors & vitamins
"""
import csv
from pathlib import Path

PATHWAYS = {
    "GLY": "Glycolysis", "PPP": "Pentose phosphate pathway", "TCA": "TCA cycle", "GLN": "Glutamine metabolism",
    "GLU": "Glutamate metabolism", "AA": "Amino acid metabolism", "FA": "Fatty acid metabolism",
    "PUR": "Purine metabolism", "PYR": "Pyrimidine metabolism", "OC": "One-carbon metabolism",
    "RED": "Redox metabolism", "NUC": "Nucleotide biosynthesis", "HEX": "Hexosamine & nucleotide sugars",
    "PL": "Phospholipid head groups", "UREA": "Urea cycle & polyamines", "COF": "Cofactors & vitamins",
}

# id, standard name, label, synonyms (|), formula, KEGG, HMDB, PubChem CID, ChEBI, pathways (|), note
L = [
 # --- glycolysis / gluconeogenesis
 ("glc", "D-Glucose", "Glucose", "Glucose|D-(+)-Glucose|Dextrose|alpha-D-Glucose|beta-D-Glucose", "C6H12O6", "C00031", "HMDB0000122", "5793", "4167", "GLY", "Hexose isomers (galactose, mannose, fructose) share the formula; confirm by RT."),
 ("g6p", "Glucose 6-phosphate", "G6P", "D-Glucose 6-phosphate|Glucose-6-phosphate|G6P", "C6H13O9P", "C00092", "HMDB0001401", "5958", "4170", "GLY|PPP|HEX", "Hexose-phosphate isomers (F6P, M6P, G1P) share the formula."),
 ("f6p", "Fructose 6-phosphate", "F6P", "D-Fructose 6-phosphate|Fructose-6-phosphate|F6P", "C6H13O9P", "C00085", "HMDB0000124", "69507", "", "GLY|HEX", ""),
 ("fbp", "Fructose 1,6-bisphosphate", "FBP", "D-Fructose 1,6-bisphosphate|Fructose-1,6-bisphosphate|FBP", "C6H14O12P2", "C00354", "HMDB0001058", "718", "", "GLY", ""),
 ("dhap", "Dihydroxyacetone phosphate", "DHAP", "Glycerone phosphate|DHAP", "C3H7O6P", "C00111", "HMDB0001473", "668", "", "GLY", ""),
 ("ga3p", "Glyceraldehyde 3-phosphate", "GA3P", "D-Glyceraldehyde 3-phosphate|GAP|G3P", "C3H7O6P", "C00118", "HMDB0001112", "729", "", "GLY|PPP", "Same formula as DHAP."),
 ("3pg", "3-Phosphoglycerate", "3PG", "3-Phosphoglyceric acid|3-Phospho-D-glycerate|3PG", "C3H7O7P", "C00197", "HMDB0000807", "724", "", "GLY|OC", "Same formula as 2PG."),
 ("2pg", "2-Phosphoglycerate", "2PG", "2-Phosphoglyceric acid|2-Phospho-D-glycerate|2PG", "C3H7O7P", "C00631", "HMDB0000362", "59", "", "GLY", "Same formula as 3PG; often co-elute."),
 ("pep", "Phosphoenolpyruvate", "PEP", "Phosphoenolpyruvic acid|PEP", "C3H5O6P", "C00074", "HMDB0000263", "1005", "", "GLY", ""),
 ("pyr", "Pyruvate", "Pyruvate", "Pyruvic acid", "C3H4O3", "C00022", "HMDB0000243", "1060", "15361", "GLY|TCA", ""),
 ("lac", "L-Lactate", "Lactate", "Lactate|Lactic acid|L-Lactic acid|DL-Lactic acid", "C3H6O3", "C00186", "HMDB0000190", "107689", "422", "GLY", ""),
 ("gro3p", "Glycerol 3-phosphate", "Glycerol-3P", "sn-Glycerol 3-phosphate|Glycerol-3-phosphate|Glycerophosphoric acid", "C3H9O6P", "C00093", "HMDB0000126", "439162", "15978", "GLY|PL", ""),
 ("g1p", "Glucose 1-phosphate", "G1P", "alpha-D-Glucose 1-phosphate|α-Glucose-1-phosphate|Glucose-1-phosphate", "C6H13O9P", "C00103", "HMDB0001586", "439165", "", "GLY|HEX", ""),
 ("g16bp", "Glucose 1,6-bisphosphate", "G1,6BP", "α-D-Glucose-1,6-bisphosphate|alpha-D-Glucose 1,6-bisphosphate", "C6H14O12P2", "C01231", "", "", "", "GLY", "Same formula as FBP."),
 ("m6p", "Mannose 6-phosphate", "M6P", "D-Mannose 6-phosphate|Mannose-6-phosphate", "C6H13O9P", "C00275", "HMDB0001078", "65127", "", "HEX", "Hexose-phosphate isomer; separation from G6P/F6P not guaranteed."),
 ("glyc", "Glycerate", "Glycerate", "Glyceric acid|D-Glyceric acid|D-Glycerate", "C3H6O4", "C00258", "HMDB0000139", "439194", "", "GLY|OC", ""),
 # --- pentose phosphate pathway
 ("6pg", "6-Phosphogluconate", "6PG", "6-Phospho-D-gluconate|6-Phosphogluconic acid", "C6H13O10P", "C00345", "HMDB0001316", "91493", "", "PPP", ""),
 ("glcn", "Gluconate", "Gluconate", "Gluconic acid|D-Gluconic acid|D-Gluconate", "C6H12O7", "C00257", "HMDB0000625", "10690", "", "PPP", "Same formula as galactonate."),
 ("ru5p", "Ribulose 5-phosphate", "Ru5P", "D-Ribulose 5-phosphate|Ribulose-5-phosphate", "C5H11O8P", "C00199", "HMDB0000618", "439184", "", "PPP", "Pentose-phosphate isomers (R5P, Xu5P) share the formula."),
 ("r5p", "Ribose 5-phosphate", "R5P", "D-Ribose 5-phosphate|Ribose-5-phosphate", "C5H11O8P", "C00117", "HMDB0001548", "", "", "PPP|NUC", ""),
 ("xu5p", "Xylulose 5-phosphate", "Xu5P", "D-Xylulose 5-phosphate", "C5H11O8P", "C00231", "HMDB0000868", "439190", "", "PPP", ""),
 ("s7p", "Sedoheptulose 7-phosphate", "S7P", "D-Sedoheptulose 7-phosphate|Sedoheptulose-7-phosphate", "C7H15O10P", "C05382", "HMDB0001068", "165007", "", "PPP", ""),
 ("e4p", "Erythrose 4-phosphate", "E4P", "D-Erythrose 4-phosphate", "C4H9O7P", "C00279", "HMDB0001321", "122357", "", "PPP", ""),
 ("prpp", "PRPP", "PRPP", "5-Phosphoribosyl 1-pyrophosphate|Phosphoribosyl pyrophosphate", "C5H13O14P3", "C00119", "HMDB0000280", "7339", "", "PPP|NUC", ""),
 # --- TCA cycle and related
 ("cit", "Citrate", "Citrate", "Citric acid", "C6H8O7", "C00158", "HMDB0000094", "311", "30769", "TCA", "Isocitrate has the same formula; usually separated chromatographically."),
 ("acon", "cis-Aconitate", "Aconitate", "Aconitic acid|cis-Aconitic acid|trans-Aconitic acid", "C6H6O6", "C00417", "HMDB0000072", "643757", "", "TCA", ""),
 ("akg", "2-Oxoglutarate", "α-KG", "alpha-Ketoglutarate|alpha-Ketoglutaric acid|2-Oxoglutaric acid|α-Ketoglutaric acid|2-Ketoglutaric acid|Oxoglutaric acid", "C5H6O5", "C00026", "HMDB0000208", "51", "", "TCA|GLU", ""),
 ("suc", "Succinate", "Succinate", "Succinic acid", "C4H6O4", "C00042", "HMDB0000254", "1110", "15741", "TCA", "Same formula as methylmalonate."),
 ("fum", "Fumarate", "Fumarate", "Fumaric acid", "C4H4O4", "C00122", "HMDB0000134", "444972", "", "TCA", "Maleate (cis isomer) has the same formula."),
 ("mal", "Malate", "Malate", "Malic acid|L-Malic acid|L-(-)-Malic acid|(S)-Malate", "C4H6O5", "C00149", "HMDB0000156", "222656", "", "TCA", ""),
 ("oaa", "Oxaloacetate", "OAA", "Oxaloacetic acid", "C4H4O5", "C00036", "HMDB0000223", "970", "", "TCA", "Unstable; rarely measured directly."),
 ("ita", "Itaconate", "Itaconate", "Itaconic acid", "C5H6O4", "C00490", "HMDB0002092", "811", "", "TCA", "Made from cis-aconitate by ACOD1/IRG1."),
 # --- glutamine / glutamate
 ("gln", "L-Glutamine", "Glutamine", "Glutamine", "C5H10N2O3", "C00064", "HMDB0000641", "5961", "18050", "GLN|AA", ""),
 ("glu", "L-Glutamate", "Glutamate", "Glutamate|Glutamic acid|L-Glutamic acid", "C5H9NO4", "C00025", "HMDB0000148", "33032", "16015", "GLU|GLN|AA", ""),
 ("pglu", "5-Oxoproline", "Pyroglutamate", "Pyroglutamic acid|L-Pyroglutamic acid|5-Oxo-L-proline|Pidolic acid", "C5H7NO3", "C01879", "HMDB0000267", "7405", "", "GLU|RED", "Also forms in-source from glutamate/glutamine."),
 ("gaba", "4-Aminobutanoate", "GABA", "gamma-Aminobutyric acid|GABA|4-Aminobutyric acid", "C4H9NO2", "C00334", "HMDB0000112", "119", "16865", "GLU", ""),
 ("naa", "N-Acetyl-L-aspartate", "NAA", "N-Acetylaspartic acid|N-Acetylaspartate|NAA", "C6H9NO5", "C01042", "HMDB0000812", "65065", "", "AA|GLU", ""),
 ("naag", "N-Acetylaspartylglutamate", "NAAG", "N-Acetyl-1-aspartylglutamic acid|NAAG", "C11H16N2O8", "C12270", "HMDB0001067", "5255", "", "GLU", ""),
 # --- amino acids
 ("asp", "L-Aspartate", "Aspartate", "Aspartic acid|L-Aspartic acid|Aspartate", "C4H7NO4", "C00049", "HMDB0000191", "5960", "17053", "AA|NUC|TCA", ""),
 ("asn", "L-Asparagine", "Asparagine", "Asparagine", "C4H8N2O3", "C00152", "HMDB0000168", "6267", "17196", "AA", ""),
 ("asp4p", "4-Phospho-L-aspartate", "Asp-4P", "L-Aspartyl-4-phosphate|Aspartyl phosphate", "C4H8NO7P", "C03082", "", "", "", "AA", "Not a typical mammalian metabolite; check annotation."),
 ("ser", "L-Serine", "Serine", "Serine", "C3H7NO3", "C00065", "HMDB0000187", "5951", "17115", "AA|OC", ""),
 ("gly", "Glycine", "Glycine", "", "C2H5NO2", "C00037", "HMDB0000123", "750", "15428", "AA|OC", ""),
 ("ala", "L-Alanine", "Alanine", "Alanine", "C3H7NO2", "C00041", "HMDB0000161", "5950", "16977", "AA", ""),
 ("pro", "L-Proline", "Proline", "Proline", "C5H9NO2", "C00148", "HMDB0000162", "145742", "17203", "AA", ""),
 ("phe", "L-Phenylalanine", "Phe", "Phenylalanine", "C9H11NO2", "C00079", "HMDB0000159", "6140", "17295", "AA", "Essential amino acid: not made from glucose."),
 ("tyr", "L-Tyrosine", "Tyr", "Tyrosine", "C9H11NO3", "C00082", "HMDB0000158", "6057", "17895", "AA", ""),
 ("trp", "L-Tryptophan", "Trp", "Tryptophan", "C11H12N2O2", "C00078", "HMDB0000929", "6305", "16828", "AA", "Essential amino acid."),
 ("ile", "L-Isoleucine", "Ile", "Isoleucine", "C6H13NO2", "C00407", "HMDB0000172", "6306", "17191", "AA", "Leucine/alloisoleucine share the formula."),
 ("leu", "L-Leucine", "Leu", "Leucine", "C6H13NO2", "C00123", "HMDB0000687", "6106", "15603", "AA", ""),
 ("val", "L-Valine", "Val", "Valine", "C5H11NO2", "C00183", "HMDB0000883", "6287", "16414", "AA", ""),
 ("lys", "L-Lysine", "Lys", "Lysine", "C6H14N2O2", "C00047", "HMDB0000182", "5962", "18019", "AA", ""),
 ("his", "L-Histidine", "His", "Histidine", "C6H9N3O2", "C00135", "HMDB0000177", "6274", "15971", "AA", ""),
 ("met", "L-Methionine", "Methionine", "Methionine|L-(-)-Methionine", "C5H11NO2S", "C00073", "HMDB0000696", "6137", "16643", "AA|OC", "Essential amino acid."),
 ("mhis1", "1-Methylhistidine", "1-MeHis", "", "C7H11N3O2", "C01152", "HMDB0000001", "92105", "", "AA", ""),
 ("uro", "Urocanate", "Urocanate", "Urocanic acid|trans-Urocanic acid", "C6H6N2O2", "C00785", "HMDB0000301", "736715", "", "AA", ""),
 ("kyn", "L-Kynurenine", "Kynurenine", "Kynurenine", "C10H12N2O3", "C00328", "HMDB0000684", "161166", "16946", "AA|RED", "Tryptophan → NAD+ route."),
 ("aad", "2-Aminoadipate", "2-AAA", "2-Aminoadipic acid|L-2-Aminoadipic acid|alpha-Aminoadipic acid", "C6H11NO4", "C00956", "HMDB0000510", "469", "", "AA", "Lysine degradation."),
 ("aclys", "N6-Acetyl-L-lysine", "N6-AcLys", "N6-Acetyllysine", "C8H16N2O3", "C02727", "HMDB0000206", "92832", "", "AA", ""),
 ("hser", "L-Homoserine", "Homoserine", "Homoserine", "C4H9NO3", "C00263", "HMDB0000719", "12647", "", "AA", "Same formula as threonine."),
 ("glutar", "Glutarate", "Glutarate", "Glutaric acid", "C5H8O4", "C00489", "HMDB0000661", "743", "", "AA|FA", "Lysine/tryptophan degradation."),
 ("mma", "Methylmalonate", "MMA", "Methylmalonic acid", "C4H6O4", "C02170", "HMDB0000202", "487", "", "AA|FA", "Same formula as succinate."),
 ("hcit", "Homocitrulline", "Homocitrulline", "", "C7H15N3O3", "C02427", "HMDB0000679", "65072", "", "UREA", ""),
 ("hip", "Hippurate", "Hippurate", "Hippuric acid|Benzoylglycine", "C9H9NO3", "C01586", "HMDB0000714", "464", "18089", "AA", "Benzoate conjugate; largely diet/microbiome-derived in vivo."),
 ("tau", "Taurine", "Taurine", "", "C2H7NO3S", "C00245", "HMDB0000251", "1123", "15891", "AA|RED", ""),
 ("htau", "Hypotaurine", "Hypotaurine", "", "C2H7NO2S", "C00519", "HMDB0000965", "107812", "16668", "AA|RED", ""),
 ("crn", "Creatine", "Creatine", "", "C4H9N3O2", "C00300", "HMDB0000064", "586", "16919", "AA|UREA", ""),
 ("crtn", "Creatinine", "Creatinine", "", "C4H7N3O", "C00791", "HMDB0000562", "588", "16737", "AA", ""),
 # --- urea cycle, polyamines
 ("arg", "L-Arginine", "Arginine", "Arginine|L-(+)-Arginine", "C6H14N4O2", "C00062", "HMDB0000517", "6322", "16467", "UREA|AA", ""),
 ("orn", "L-Ornithine", "Ornithine", "Ornithine|L(+)-Ornithine", "C5H12N2O2", "C00077", "HMDB0000214", "6262", "15729", "UREA|AA", ""),
 ("citr", "L-Citrulline", "Citrulline", "Citrulline|L-(+)-Citrulline", "C6H13N3O3", "C00327", "HMDB0000904", "9750", "16349", "UREA", ""),
 ("argsuc", "Argininosuccinate", "Arg-succinate", "Argininosuccinic acid|L-Argininosuccinate", "C10H18N4O6", "C03406", "HMDB0000052", "16950", "", "UREA", ""),
 ("spd", "Spermidine", "Spermidine", "", "C7H19N3", "C00315", "HMDB0001257", "1102", "16610", "UREA", ""),
 ("acspd", "N1-Acetylspermidine", "N1-AcSpd", "", "C9H21N3O", "C00612", "HMDB0001276", "496", "", "UREA", ""),
 # --- one-carbon / methionine cycle / transsulfuration
 ("sam", "S-Adenosyl-L-methionine", "SAM", "S-Adenosylmethionine|AdoMet|SAMe", "C15H22N6O5S", "C00019", "HMDB0001185", "34755", "15414", "OC", ""),
 ("sah", "S-Adenosyl-L-homocysteine", "SAH", "S-Adenosylhomocysteine|AdoHcy", "C14H20N6O5S", "C00021", "HMDB0000939", "439155", "16680", "OC", ""),
 ("cysta", "L-Cystathionine", "Cystathionine", "Cystathionine", "C7H14N2O4S", "C02291", "HMDB0000099", "439258", "", "OC|RED|AA", "Transsulfuration."),
 ("metso", "L-Methionine sulfoxide", "Met-SO", "Methionine sulfoxide", "C5H11NO3S", "C02989", "HMDB0002005", "", "", "RED", "Same formula as 2-hydroxy-methionine."),
 # --- redox
 ("gsh", "Glutathione", "GSH", "L-Glutathione (reduced)|Glutathione (reduced)|Reduced glutathione|GSH", "C10H17N3O6S", "C00051", "HMDB0000125", "124886", "16856", "RED|GLU", ""),
 ("gssg", "Glutathione disulfide", "GSSG", "Oxidized glutathione|L-Glutathione oxidized|GSSG", "C20H32N6O12S2", "C00127", "HMDB0003337", "65359", "17858", "RED", ""),
 ("asc", "L-Ascorbate", "Ascorbate", "Ascorbic acid|Vitamin C", "C6H8O6", "C00072", "HMDB0000044", "54670067", "", "RED|COF", ""),
 ("nad", "NAD+", "NAD⁺", "Nicotinamide adenine dinucleotide (NAD+)|NAD|Nicotinamide adenine dinucleotide", "C21H27N7O14P2", "C00003", "HMDB0000902", "5893", "15846", "RED|NUC|COF", ""),
 ("nam", "Nicotinamide", "Nicotinamide", "Niacinamide", "C6H6N2O", "C00153", "HMDB0001406", "936", "17154", "RED|COF", ""),
 ("cadpr", "Cyclic ADP-ribose", "cADPR", "cADPR", "C15H21N5O13P2", "C13050", "", "", "", "RED|NUC", "NAD+ consumption product (CD38)."),
 # --- purines
 ("ade", "Adenine", "Adenine", "", "C5H5N5", "C00147", "HMDB0000034", "190", "16708", "PUR", ""),
 ("ado", "Adenosine", "Adenosine", "", "C10H13N5O4", "C00212", "HMDB0000050", "60961", "16335", "PUR", ""),
 ("amp", "AMP", "AMP", "Adenosine 5'-monophosphate|Adenosine monophosphate|Adenylic acid", "C10H14N5O7P", "C00020", "HMDB0000045", "6083", "16027", "PUR|NUC", ""),
 ("adp", "ADP", "ADP", "Adenosine diphosphate (ADP)|Adenosine 5'-diphosphate|Adenosine diphosphate", "C10H15N5O10P2", "C00008", "HMDB0001341", "6022", "16761", "PUR|NUC", ""),
 ("atp", "ATP", "ATP", "Adenosine 5'-triphosphate|Adenosine triphosphate", "C10H16N5O13P3", "C00002", "HMDB0000538", "5957", "15422", "PUR|NUC", ""),
 ("imp", "IMP", "IMP", "Inosinic acid|Inosine-5'-monophosphate (IMP)|Inosine monophosphate|Inosine 5'-monophosphate", "C10H13N4O8P", "C00130", "HMDB0000175", "135398640", "17202", "PUR|NUC", ""),
 ("ino", "Inosine", "Inosine", "", "C10H12N4O5", "C00294", "HMDB0000195", "135398641", "17596", "PUR", "Allopurinol riboside has the same formula."),
 ("hyp", "Hypoxanthine", "Hypoxanthine", "", "C5H4N4O", "C00262", "HMDB0000157", "135398638", "17368", "PUR", ""),
 ("xan", "Xanthine", "Xanthine", "", "C5H4N4O2", "C00385", "HMDB0000292", "1188", "17712", "PUR", ""),
 ("urate", "Urate", "Urate", "Uric acid", "C5H4N4O3", "C00366", "HMDB0000289", "1175", "17775", "PUR|RED", ""),
 ("gua", "Guanine", "Guanine", "", "C5H5N5O", "C00242", "HMDB0000132", "135398634", "16235", "PUR", "Same formula as 2-hydroxyadenine."),
 ("s-amp", "Adenylosuccinate", "S-AMP", "N(6)-(1,2-dicarboxyethyl)-AMP|Adenylosuccinic acid|Succinyl-AMP", "C14H18N5O11P", "C03794", "HMDB0000536", "", "", "PUR|NUC", ""),
 ("sado", "Succinyladenosine", "S-Ado", "N6-Succinyl Adenosine|Succinyladenosine", "C14H17N5O8", "", "HMDB0000912", "", "", "PUR", "Dephosphorylated adenylosuccinate."),
 # --- pyrimidines
 ("oro", "Orotate", "Orotate", "Orotic acid", "C5H4N2O4", "C00295", "HMDB0000226", "967", "16742", "PYR|NUC", ""),
 ("orod", "Orotidine", "Orotidine", "", "C10H12N2O8", "", "HMDB0000788", "", "", "PYR|NUC", "Nucleoside of OMP."),
 ("ura", "Uracil", "Uracil", "", "C4H4N2O2", "C00106", "HMDB0000300", "1174", "17568", "PYR", ""),
 ("urd", "Uridine", "Uridine", "", "C9H12N2O6", "C00299", "HMDB0000296", "6029", "16704", "PYR", ""),
 ("ump", "UMP", "UMP", "Uridine 5'-monophosphate|Uridine monophosphate|Uridylic acid|Uridine 3'-monophosphate", "C9H13N2O9P", "C00105", "HMDB0000288", "6030", "16695", "PYR|NUC", "3'- and 5'-UMP share the formula; the cellular pool is mainly 5'-UMP."),
 ("udp", "UDP", "UDP", "Uridine 5'-diphosphate|Uridine diphosphate", "C9H14N2O12P2", "C00015", "HMDB0000295", "6031", "17659", "PYR|NUC", ""),
 ("cmp", "CMP", "CMP", "Cytidine 5'-monophosphate|Cytidine 5'-monophosphate (hydrate)|Cytidine monophosphate|Cytidylic acid", "C9H14N3O8P", "C00055", "HMDB0000095", "6131", "17361", "PYR|NUC", ""),
 ("cyt", "Cytosine", "Cytosine", "", "C4H5N3O", "C00380", "HMDB0000630", "597", "16040", "PYR", ""),
 # --- hexosamine & nucleotide sugars
 ("glcn6p", "Glucosamine 6-phosphate", "GlcN6P", "D-Glucosamine 6-phosphate|Glucosamine-6-phosphate", "C6H14NO8P", "C00352", "HMDB0001254", "439217", "", "HEX", ""),
 ("udpglcnac", "UDP-N-acetylglucosamine", "UDP-GlcNAc", "UDP-GlcNAc|Uridine diphosphate N-acetylglucosamine", "C17H27N3O17P2", "C00043", "HMDB0000290", "445675", "16264", "HEX|PYR", "UDP-GalNAc has the same formula."),
 ("udpgal", "UDP-galactose", "UDP-Gal", "Uridine 5'-diphosphogalactose|UDP-D-galactose", "C15H24N2O17P2", "C00052", "HMDB0000302", "18068", "", "HEX|PYR", "UDP-glucose has the same formula; usually co-elute."),
 ("neu5ac", "N-Acetylneuraminate", "Neu5Ac", "N-Acetylneuraminic acid|Sialic acid|Neu5Ac", "C11H19NO9", "C19910", "HMDB0000230", "445063", "", "HEX", ""),
 ("neu5gc", "N-Glycolylneuraminate", "Neu5Gc", "N-Glycolylneuraminic acid|Neu5Gc", "C11H19NO10", "C03410", "HMDB0000833", "", "", "HEX", "Not made by human cells; usually taken up from serum."),
 # --- phospholipid head groups
 ("cho", "Choline", "Choline", "", "C5H13NO", "C00114", "HMDB0000097", "305", "15354", "PL|OC", ""),
 ("gpc", "Glycerophosphocholine", "GPC", "sn-Glycero-3-phosphocholine|Glycerylphosphorylcholine|alpha-GPC", "C8H20NO6P", "C00670", "HMDB0000086", "71920", "", "PL", ""),
 ("pe", "Phosphoethanolamine", "PEtn", "O-Phosphoethanolamine|Ethanolamine phosphate", "C2H8NO4P", "C00346", "HMDB0000224", "1015", "", "PL", ""),
 ("gpe", "Glycerophosphoethanolamine", "GPE", "Glycerylphosphorylethanolamine|sn-Glycero-3-phosphoethanolamine", "C5H14NO6P", "C01233", "HMDB0000114", "", "", "PL", ""),
 ("cdpcho", "CDP-choline", "CDP-Cho", "Cytidine 5'-diphosphocholine|Citicoline|CDP-Choline", "C14H26N4O11P2", "C00307", "HMDB0001413", "13804", "", "PL|PYR", ""),
 ("cdpetn", "CDP-ethanolamine", "CDP-Etn", "CDP-Ethanolamine", "C11H20N4O11P2", "C00570", "HMDB0001564", "", "", "PL|PYR", ""),
 ("gpg", "Glycerophosphoglycerol", "GPG", "", "C6H15O8P", "C03274", "", "", "", "PL", ""),
 ("gpi", "Glycerophosphoinositol", "GPI", "sn-Glycero-3-phospho-1-inositol", "C9H19O11P", "C01225", "HMDB0011649", "", "", "PL", ""),
 ("ach", "Acetylcholine", "ACh", "", "C7H15NO2", "C01996", "HMDB0000895", "187", "15355", "PL", ""),
 # --- fatty acid metabolism / carnitine
 ("car", "L-Carnitine", "Carnitine", "Carnitine", "C7H15NO3", "C00318", "HMDB0000062", "10917", "16347", "FA", ""),
 ("c2car", "Acetylcarnitine", "C2-carnitine", "Acetyl-L-carnitine|O-Acetylcarnitine|ALCAR", "C9H17NO4", "C02571", "HMDB0000201", "7045767", "", "FA", "Acetyl group can carry glucose-derived carbon (M+2)."),
 ("c4car", "Butyrylcarnitine", "C4-carnitine", "", "C11H21NO4", "C02862", "HMDB0002013", "", "", "FA", ""),
 ("c8car", "Octanoylcarnitine", "C8-carnitine", "L-Octanoylcarnitine", "C15H29NO4", "C02838", "HMDB0000791", "", "", "FA", ""),
 ("c10car", "Decanoylcarnitine", "C10-carnitine", "", "C17H33NO4", "C03299", "HMDB0000651", "", "", "FA", ""),
 ("c16car", "Palmitoylcarnitine", "C16-carnitine", "L-Palmitoylcarnitine", "C23H45NO4", "C02990", "HMDB0000222", "", "", "FA", ""),
 ("bhb", "3-Hydroxybutyrate", "3-HB", "3-Hydroxybutyric acid|beta-Hydroxybutyrate|(R)-3-Hydroxybutanoate", "C4H8O3", "C01089", "HMDB0000011", "441", "", "FA", "Ketone body; same formula as 4-hydroxybutyrate."),
 ("ala_ln", "alpha-Linolenic acid", "ALA", "α-Linolenic acid|Linolenic acid", "C18H30O2", "C06427", "HMDB0001388", "5280934", "", "FA", "Essential fatty acid."),
 ("malon", "Malonate", "Malonate", "Malonic acid", "C3H4O4", "C00383", "HMDB0000691", "867", "", "FA", ""),
 ("adip", "Adipate", "Adipate", "Adipic acid", "C6H10O4", "C06104", "HMDB0000448", "196", "", "FA", "Dicarboxylic acid (ω-oxidation)."),
 ("suber", "Suberate", "Suberate", "Suberic acid", "C8H14O4", "C08278", "HMDB0000893", "10457", "", "FA", ""),
 # --- cofactors / vitamins and other
 ("pant", "Pantothenate", "Pantothenate", "Pantothenic acid|Vitamin B5", "C9H17NO5", "C00864", "HMDB0000210", "6613", "", "COF", "Vitamin; not synthesized by human cells."),
 ("threon", "Threonate", "Threonate", "Threonic acid|L-Threonic acid", "C4H8O5", "C01620", "HMDB0000943", "", "", "RED", "Ascorbate breakdown product."),
 ("glucar", "Glucarate", "Glucarate", "Glucaric acid|D-Glucaric acid|Saccharic acid", "C6H10O8", "C00818", "HMDB0000663", "", "", "HEX", ""),
 ("galacton", "Galactonate", "Galactonate", "Galactonic acid|D-Galactonate", "C6H12O7", "C00880", "HMDB0000565", "", "", "HEX", "Same formula as gluconate."),
 ("galur", "Galacturonate", "GalUA", "Galacturonic acid|D-Galacturonic acid", "C6H10O7", "C00333", "HMDB0002545", "", "", "HEX", "Same formula as glucuronate."),
 ("arabitol", "Arabitol", "Arabitol", "D-Arabitol", "C5H12O5", "C01904", "HMDB0001851", "", "", "PPP", "Pentitol; ribitol/xylitol share the formula."),
]

# probable mis-annotations from search libraries: data name -> (library id, reason). Shown as LOW confidence,
# never assigned automatically; the user can accept them in the annotation editor.
SUGGESTIONS = {
    "3-Oxoglutaric acid": ("akg", "Same formula as 2-oxoglutarate (α-KG); 3-oxoglutarate is not a common cellular metabolite."),
    "Maleic acid": ("fum", "Maleate is the cis isomer of fumarate; in cells the TCA intermediate is fumarate."),
    "4-Oxoproline": ("pglu", "Same formula as 5-oxoproline (pyroglutamate)."),
    "Glutaric anhydride": ("glutar", "Anhydride is a chemical artifact; likely glutarate-related in-source species."),
    "2-Hydroxyadenine": ("gua", "Same formula as guanine."),
    "Allopurinol riboside": ("ino", "Same formula as inosine (drug metabolite unlikely in culture)."),
    "D-myo-Inositol 1,4-bisphosphate": ("fbp", "Same formula as hexose bisphosphates."),
    "Heptulose-2-Phosphate": ("s7p", "Same formula as sedoheptulose 7-phosphate."),
    "4-Hydroxybutyric acid (GHB)": ("bhb", "Same formula as 3-hydroxybutyrate."),
    "2-Hydroxy-L-methionine": ("metso", "Same formula as methionine sulfoxide."),
    "DL-β-Leucine": ("leu", "Same formula as leucine/isoleucine."),
    "L-Alloisoleucine": ("ile", "Same formula as isoleucine/leucine."),
}

out = Path(__file__).with_name("metabolite_library.csv")
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["LibraryID", "StandardName", "Label", "Synonyms", "Formula", "KEGG", "HMDB", "PubChemCID", "ChEBI",
                "Pathways", "Note", "IDSource"])
    for r in L:
        w.writerow(list(r) + ["curated (not machine-verified)"])
sg = Path(__file__).with_name("annotation_suggestions.csv")
with open(sg, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["DataName", "LibraryID", "Reason"])
    for k, (i, why) in SUGGESTIONS.items():
        w.writerow([k, i, why])
pw = Path(__file__).with_name("pathways.csv")
with open(pw, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["Code", "Pathway"])
    for k, v in PATHWAYS.items():
        w.writerow([k, v])
print(len(L), "library entries")
