#!/usr/bin/env python3

import sys
import cgi
import re, os, shutil
from collections import OrderedDict
import subprocess
import time


##### v1.1.7 update
# Presence or abusence of databases are checked only on mode E.
# The function check_presense_of_databases was modified.

# The BHnum_Xxxxx-xxx columns in the result.csv index line are now made from taxon sapling list.
# The function make_lines_atmarkSeparated was modified.
######

##### v1.1.6 update
# Lengths of species names are imposed a limitation. It should be less than or equal to 38.
# The function check_pickup_taxonSampling was modified.
######

AddintHeaderAfterAT = "D"   ## L:leave or D:Delete @xxxx for the summarize analysis.
draw_speciesTree = "Not"  ## Draw: Draw species tree in the .pdf file. or Not:

if len (sys.argv) < 2:
    print("Error. you need two arguments.")
    print("Example: ./main_script.py ENSORLT00000007282.1")
    exit()

queryID = sys.argv[1]

###################
geneticCode = {
            "CTA"                             : "L",
            "CTT"                             : "L",
            "CTG"                             : "L",
            "CTC"                             : "L",
            "TTA"                             : "L",
            "TTG"                             : "L",
            "TTR"                             : "L",
           # "((CG.)|(AG(A|G|R)))"               : "R",
            "CGA"                             : "R",
            "CGT"                             : "R",
            "CGG"                             : "R",
            "CGC"                             : "R",
            "AGA"                             : "R",
            "AGG"                             : "R",
            "AGR"                             : "R",
           # "(((U|T)C.)|(AG(U|T|C|Y)))"         : "S",
            "TCA"                             : "S",
            "TCT"                             : "S",
            "TCG"                             : "S",
            "TCC"                             : "S",
            "AGT"                             : "S",
            "AGC"                             : "S",
            "AGY"                             : "S",
           # "(GC.)"                             : "A",
            "GCA"                             : "A",
            "GCT"                             : "A",
            "GCG"                             : "A",
            "GCC"                             : "A",
           # "(GG.)"                             : "G",
            "GGT"                             : "G",
            "GGC"                             : "G",
            "GGA"                             : "G",
            "GGG"                             : "G",
           # "(CC.)"                             : "P",
            "CCT"                             : "P",
            "CCC"                             : "P",
            "CCA"                             : "P",
            "CCG"                             : "P",
           # "(AC.)"                             : "T",
            "ACT"                             : "T",
            "ACC"                             : "T",
            "ACA"                             : "T",
            "ACG"                             : "T",
           # "(G(U|T).)"                         : "V",
            "GTT"                             : "V",
            "GTC"                             : "V",
            "GTA"                             : "V",
            "GTG"                             : "V",
           # "(A(U|T)(U|T|C|Y|A))"               : "I",
            "ATT"                             : "I",
            "ATC"                             : "I",
            "ATA"                             : "I",
           # "(((U|T)A(A|G|R))|((T|U)GA))"    : "_",
            "TAA"                             : "X",  #* Ter
            "TAG"                             : "X",  #* Ter
            "TAR"                             : "X",  #* Ter
            "TGA"                             : "X",  #* Ter
           # "((U|T)G(U|T|C|Y))"                 : "C",
            "TGT"                             : "C",  #Cys  
            "TGC"                             : "C",  #Cys  
            "TGY"                             : "C",  #Cys  
           # "(GA(U|T|C|Y))"                     : "D",
            "GAT"                             : "D",  #Asp
            "GAC"                             : "D",  #Asp
            "GAY"                             : "D",  #Asp
           # "(GA(A|G|R))"                       : "E",
            "GAA"                             : "E",  #Glu
            "GAG"                             : "E",  #Glu
            "GAR"                             : "E",  #Glu
           # "((U|T)(U|T)(U|T|C|Y))"             : "F",
            "TTT"                             : "F",  #Phe
            "TTC"                             : "F",  #Phe
            "TTY"                             : "F",  #Phe
           # "(CA(U|T|C|Y))"                     : "H",
            "CAT"                             : "H",  #His
            "CAC"                             : "H",  #His
            "CAY"                             : "H",  #His
           # "(AA(A|G|R))"                       : "K",
            "AAA"                             : "K",  #Lys
            "AAG"                             : "K",  #Lys
            "AAR"                             : "K",  #Lys
           # "(AA(U|T|C|Y))"                     : "N",
            "AAT"                             : "N",  #Asn
            "AAC"                             : "N",  #Asn
            "AAY"                             : "N",  #Asn
           # "(CA(A|G|R))"                       : "Q",
            "CAA"                             : "Q",  #Gln
            "CAG"                             : "Q",  #Gln
            "CAR"                             : "Q",  #Gln
           # "((U|T)A(U|T|C|Y))"                 : "Y",
            "TAT"                             : "Y",  #Tyr
            "TAC"                             : "Y",  #Tyr
            "TAY"                             : "Y",  #Tyr
           # "(A(U|T)G)"                         : "M",
            "ATG"                             : "M",  #Met
           # "((U|T)GG)"                         : "W",
            "TGG"                             : "W",  #Trp  
           # "..."                               : "X",
           # "(NNN)"                             : "X",
            "NNN"                             : "X",  
           # "(N(.|N).)"                         : "X",
            "N.."                             : "X",  
            "NN."                             : "X",  
           # "(.(.|N)N)"                         : "X",
            ".NN"                             : "X",  
            "..N"                             : "X",  
           # "(.N.)"                             : "X",
            ".N."                             : "X",  
            "---"                             : "-"}


resHTMLlines_2steps = '''
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>
    <head>
        <meta http-equiv="content-type" content="text/html;charset=shift_jis">
        <title>ORTHOSCOPE STAR: res 2steps</title>
        <link href="main.css" rel="stylesheet" type="text/css" media="all">
    </head>


<body bgcolor="#eeeeee" leftmargin="20" marginheight="20" marginwidth="20" topmargin="20">
<table align="center" border="0" cellspacing="5" cellpadding="5" bgcolor="white">

  <!-- title -->
  <tr><td width="600"><table width="100%" border="0" cellspacing="2" cellpadding="0" bgcolor="#000088" height="50">
     <tr>
       <td align="center" valign="middle">
         <font size="5" color=#FFFFFF face="Verdana, Arial, Helvetica, sans-serif"><b>ORTHOSCOPE STAR</b></font>
       </td>
     </tr>
  </table></td></tr>


  <!-- Result table -->
  <tr><td><table border="1">

    <tr>
      <td align="center"  valign="top" width="100">Query sequence:</td>
      <td align="center"  valign="top" width="400">FIRSTQUERY</td>
      <td align="center" valign="top" width="100">&nbsp;</td>
    </tr>
    <tr>
      <td align="center" valign="top"><!-- PREVIOUSPAGE --></td>
      <td align="center" valign="top"><a href="EACHDIRADDRESS_100_analysisSummary.txt" target="_blank">Summary</a></td>
      <td align="center" valign="top"><!-- NEXTPAGE --></td>
    </tr>

    <!--
    <tr>
      <td align="center" valign="top"><b>2nd tree</b></td>
      <td align="center" valign="top">&nbsp;</td>
      <td align="center" valign="top">&nbsp;</td>
    </tr>

    <tr>
      <td align="left" valign="top">Bs value of sister clade monophyly:</td>
      <td align="left" valign="top">&nbsp;BS_2NDTREE_SISTERNODE</td>
      <td align="right" valign="top">&nbsp;</td>
    </tr>

    <tr>
      <td align="left" valign="top">Bs value of sister vs query-gene groups:</td>
      <td align="left" valign="top">&nbsp;BS_2NDTREE_PARENTNODE</td>
      <td align="right" valign="top">&nbsp;</td>
    </tr>

    <tr>
      <td align="left" valign="top">Bs value of vertebrate clade monophyly:</td>
      <td align="left" valign="top">&nbsp;BS_2NDTREE_VERTEBRATENODE</td>
      <td align="right" valign="top">&nbsp;</td>
    </tr>

    <tr>
      <td align="center" valign="top"><b>1st tree</b></td>
      <td align="center" valign="top">&nbsp;</td>
      <td align="center" valign="top">&nbsp;</td>
    </tr>
    -->

    <tr>
      <td align="center" valign="top">Bs value of orthogroup-basal node (1st tree):</td>
      <td align="center" valign="top">&nbsp;BSVALUE_orthogroup_1STTREE</td>
      <td align="center" valign="top">&nbsp;</td>
    </tr>

    </table></td></tr>

    <tr><td><hr></td></tr>


    <!-- Tree talbe -->
    <tr><td><table width="100%" border="0">
      <tr>
        <td colspan="2"><b>2nd tree:</b> Speciation/duplication events in the query sequence lineage </td>
      </tr>
      <tr>
        <td align="left" valign="top"><img src="EACHDIRADDRESS_240_2ndRearranged_geneTree.png"></td>
        <td align="left" valign="top"><img src="EACHDIRADDRESS_240_2ndGeneTree.png"></td>
      </tr>
      <tr>
        <td align="center" valign="top">Rearranged gene tree (<a href="EACHDIRADDRESS_240_2ndRearranged_geneTree.pdf" target="_blank">PDF</a>)</td>
        <td align="center" valign="top">Gene tree (<a href="EACHDIRADDRESS_240_2ndGeneTree.pdf" target="_blank">PDF</a>)</td>
      </tr>

      <!-- 
      <tr>
        <td>Alignment: <a href="170_aln_prot.html" target="_blank">Amino acid</a>, <a href="190_aln_nucl.txt" target="_blank">Nucleotide</a></td>
        <td>&nbsp;</td>
      </tr>
      -->

      <tr><td colspan="2"><Hr></td></tr>
      <tr>
        <td colspan="2"><b>1st tree:</b> Orthogroup</td>
      </tr>
      <tr>
        <td align="left" valign="top"><img src="EACHDIRADDRESS_115_1stRearranged_geneTree.png"></td>
        <td align="left" valign="top"><img src="EACHDIRADDRESS_115_1stGeneTree.png"></td>
      </tr>

      <tr>
        <td align="center" valign="top" name="REARRANGED2">Rearranged gene tree (<a href="EACHDIRADDRESS_115_1stRearranged_geneTree.pdf" target="_blank">PDF</a>)</td>
        <td align="center" valign="top">NJ tree (<a href="EACHDIRADDRESS_115_1stGeneTree.pdf" target="_blank">PDF</a>)</td>
      </tr>

    </table></td></tr>


      <tr>
        <td colspan="2"><hr></td>
      </tr>

</table>
</body>
</html>
'''

resHTMLlines_incomplete = '''
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>
    <head>
        <meta http-equiv="content-type" content="text/html;charset=shift_jis">
        <title>ORTHOSCOPE STAR: incomplete</title>
        <link href="main.css" rel="stylesheet" type="text/css" media="all">
    </head>


<body bgcolor="#eeeeee" leftmargin="20" marginheight="20" marginwidth="20" topmargin="20">
<table align="center" border="0" cellspacing="5" cellpadding="5" bgcolor="white">

  <!-- title -->
  <tr><td width="600"><table width="100%" border="0" cellspacing="2" cellpadding="0" bgcolor="#000088" height="50">
     <tr>
       <td align="center" valign="middle">
         <font size="5" color=#FFFFFF face="Verdana, Arial, Helvetica, sans-serif"><b>ORTHOSCOPE STAR</b></font>
       </td>
     </tr>
  </table></td></tr>

  <!-- Result table -->
  <tr><td><table border="1">

    <tr>
      <td align="center"  valign="top" width="100">Query sequence:</td>
      <td align="center"  valign="top" width="400">FIRSTQUERY</td>
      <td align="center" valign="top" width="100">&nbsp;</td>
    </tr>
    <tr>
      <td align="center" valign="top"><!-- PREVIOUSPAGE --></td>
      <td align="center" valign="top"><a href="EACHDIRADDRESS_100_analysisSummary.txt" target="_blank">Summary</a></td>
      <td align="center" valign="top"><!-- NEXTPAGE --></td>
    </tr>

    <!--
    <tr>
      <td align="center" valign="top"><b>2nd tree</b></td>
      <td align="center" valign="top">&nbsp;</td>
      <td align="center" valign="top">&nbsp;</td>
    </tr>

    <tr>
      <td align="left" valign="top">Bs value of sister clade monophyly:</td>
      <td align="left" valign="top">&nbsp;BS_2NDTREE_SISTERNODE</td>
      <td align="right" valign="top">&nbsp;</td>
    </tr>

    <tr>
      <td align="left" valign="top">Bs value of sister vs query-gene groups:</td>
      <td align="left" valign="top">&nbsp;BS_2NDTREE_PARENTNODE</td>
      <td align="right" valign="top">&nbsp;</td>
    </tr>

    <tr>
      <td align="left" valign="top">Bs value of vertebrate clade monophyly:</td>
      <td align="left" valign="top">&nbsp;BS_2NDTREE_VERTEBRATENODE</td>
      <td align="right" valign="top">&nbsp;</td>
    </tr>

    <tr>
      <td align="center" valign="top"><b>1st tree</b></td>
      <td align="center" valign="top">&nbsp;</td>
      <td align="center" valign="top">&nbsp;</td>
    </tr>
    -->

    <tr>
      <td align="center" valign="top">Bs value of orthogroup-basal node (1st tree):</td>
      <td align="center" valign="top">&nbsp;BSVALUE_orthogroup_1STTREE</td>
      <td align="center" valign="top">&nbsp;</td>
    </tr>

    </table></td></tr>

    <tr><td><hr></td></tr>


    <!-- Tree talbe -->
    <tr><td><table width="100%" border="0">
      <tr>
        <td colspan="2"><b>2nd tree:</b> Speciation/duplication events in the query sequence lineage </td>
      </tr>
      <tr>
        <td align="left" valign="top"><img src="EACHDIRADDRESS_240_2ndRearranged_geneTree.png"></td>
        <td align="left" valign="top"><img src="EACHDIRADDRESS_240_2ndGeneTree.png"></td>
      </tr>
      <!-- <tr>
        <td align="center" valign="top">Rearranged gene tree (<a href="240_2ndRearranged_geneTree.pdf" target="_blank">PDF</a>)</td>
        <td align="center" valign="top">Gene tree (<a href="240_2ndGeneTree.pdf" target="_blank">PDF</a>)</td>
      </tr>
      -->

      <!-- 
      <tr>
        <td>Alignment: <a href="170_aln_prot.html" target="_blank">Amino acid</a>, <a href="190_aln_nucl.txt" target="_blank">Nucleotide</a></td>
        <td>&nbsp;</td>
      </tr>
      -->

      <tr><td colspan="2"><Hr></td></tr>
      <tr>
        <td colspan="2"><b>1st tree:</b> Orthogroup</td>
      </tr>
      <tr>
        <td align="left" valign="top"><img src="EACHDIRADDRESS_115_1stRearranged_geneTree.png"></td>
        <td align="left" valign="top"><img src="EACHDIRADDRESS_115_1stGeneTree.png"></td>
      </tr>

      <!-- <tr>
        <td align="center" valign="top" name="REARRANGED2">Rearranged gene tree (<a href="115_1stRearranged_geneTree.pdf" target="_blank">PDF</a>)</td>
        <td align="center" valign="top">NJ tree (<a href="115_1stGeneTree.pdf" target="_blank">PDF</a>)</td>
      </tr> -->

    </table></td></tr>



    <tr>
      <td colspan="2"><hr></td>
    </tr>


</table>

</body>

</html>
'''


##############################################

### File checking

def check_toolsDirectory():
    dependencies = os.listdir(path='tools')
    if mode == "E":
        if not "blastp" in dependencies:
            print("Error: Cannot find blastp in your tools directory.")
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
        if not "makeblastdb" in dependencies:
            print("Error: Cannot find makeblastdb in your tools directory.")
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
        if not "mafft" in dependencies:
            print("Error: Cannot find mafft in your tools directory.")
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
        if not "trimal" in dependencies:
            print("Error: Cannot find trimal in your tools directory.")
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
        if not "pal2nal.pl" in dependencies:
            print("Error: Cannot find pal2nal.pl in your tools directory.")
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
        if not "Notung.jar" in dependencies:
            print("Error: Cannot find Notung.jar in your tools directory.")
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
    if not "Rscript" in dependencies:
        print("Error: Cannot find Rscript in your tools directory.")
        print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
        exit()

def check_mode():
    if mode == "E" and queryID.endswith(".txt"):
        print("Error in your control.txt file.")
        print("You selected mode E. Mode E needs gene ID as the argument.")
        print("Your argument was:", queryID)
        exit()
    if mode == "D" and queryID.endswith(".txt"):
        print("Error in your control.txt file.")
        print("You selected mode D. Mode D needs gene ID as the argument.")
        print("Your argument was:", queryID)
        exit()
    if mode == "S":
        if not queryID.endswith(".txt"):
            print("Error in your control.txt file.")
            print("You selected mode S. Mode S needs .txt file as the argument.")
            print("Your argument was:", queryID)
            exit()

def check_pickup_parameter(resDict_SR, keyword):
    keyword = ">" + keyword
    parameter_SR = ""
    
    if not resDict_SR[keyword]:
       print("Error in the control file.")
       print("Check parameters in ", keyword, ".")
       exit()
    
    if not keyword in resDict_SR.keys():
        print("Cannot find ", keyword, "in your control file. Stopped.")
        exit()
    if len(resDict_SR[keyword]) > 1:
        print("Your ", keyword, " should be written in one line. Stopped.")
        exit()
    return re.sub(" ", "", resDict_SR[keyword][0])


def check_pickup_taxonSampling(dbAddress, lines_taxonSampling):
    dbLines = []
    taxonSamplingList = []
    for line in lines_taxonSampling:
        lineTMP = re.sub(" +", " ", line)
        lineTMP = re.sub("^ +", "", lineTMP)
        lineTMP = re.sub(" +$", "", lineTMP)

        line_separated = lineTMP.split(" ")
        if len (line_separated) != 3:
            print("Error in your TaxonSampling in the control file.")
            print("Each line should consists of three part separated by spaces.")
            print("Check the following line:")
            print(re.sub(" +", " ", line))
            exit()
        speciesName_color = line_separated[0]
        name_protDB_file = line_separated[1]
        name_nuclDB_file = line_separated[2]
        
        speciesName_color_separated = speciesName_color.split("_")
        if len (speciesName_color_separated) != 2:
            print("Error in your TaxonSampling in the control file.")
            print('Each species name is followed by the assigned color separated by the underscore, "_".')
            print("Check the following species_color part:")
            print(speciesName_color)
            exit()

        speciesName = speciesName_color_separated[0]
        #print("speciesName", speciesName)
        #print("len(speciesName)", len(speciesName))
        if len(speciesName) > 38:
            print("Error in your TaxonSampling in the control file.")
            print('Each species name should be less than 39 characters.')
            print("Check the following species name:")
            print(speciesName)
            exit()

        color  = speciesName_color_separated[1]
        if not color in ["Green", "Purple", "Orange", "Magenta", "Blue", "Red", "Black"]:
            print("Error in your assinged color in the control file:")
            print(speciesName_color)
            print("Assigned colors should be : Black, Green, Purple, Orange, Magenta, Blue, or Red.")
            exit()

        dbLines.append([speciesName + "_", name_protDB_file, name_nuclDB_file])
        taxonSamplingList.append(speciesName_color)

    return dbLines, taxonSamplingList


def check_presense_of_databases():
    for dbLine in dbLines:
        #print("dbLine", dbLine)
        name_protDB_file = dbLine[1]
        name_nuclDB_file = dbLine[2]
        if not os.path.isfile(dbAddress + name_protDB_file):
            print("Error in your directory:", dbAddress[:-1])
            print(name_protDB_file)
            print("is not found.")
            exit()
        if not os.path.isfile(dbAddress + name_nuclDB_file):
            print("Error in your directory:", dbAddress[:-1])
            print(name_nuclDB_file)
            print("is not found.")
            exit()


def check_controlFile(resDict_SR):
    nameLines = [">QuerySpecies",">Mode",">TaxonSampling",">SpeciesTree",">Num_rootSequences",">KeyNode",">BLAST_Evalue_threshold_for_reported_sequences",">Number_of_hits_to_report_per_genome",">ShortSequence_threshold",">Dataset",">BSthreshold",">Outdir",">Database"]
    for nameLine in nameLines:
        if not nameLine in resDict_SR.keys():
            print("Error in your control file. Nameline", nameLine, "is not found.")
            exit()


def read_controlFile():
    resDict_SR = readRes_dict("control.txt")        
    check_controlFile(resDict_SR)

    ## the others
    SpeciesTree = check_pickup_parameter(resDict_SR, "SpeciesTree")
    blastEvalue = check_pickup_parameter(resDict_SR, "BLAST_Evalue_threshold_for_reported_sequences")
    if re.search("[^\de-]", blastEvalue):
        print("Error. >BLAST_Evalue_threshold_for_reported_sequences should be floating point expression (e.g., 1e-3) or just 1.")
        print("Your >BLAST_Evalue_threshold_for_reported_sequences is", blastEvalue)
        exit()

    Number_of_hits_to_report_per_genome = check_pickup_parameter(resDict_SR, "Number_of_hits_to_report_per_genome")
    if re.search("[^\d]", Number_of_hits_to_report_per_genome):
        print("Error. >Number_of_hits_to_report_per_genome should be integer.")
        print("Your >Num_rootSequences is", Number_of_hits_to_report_per_genome)
        exit()

    Aligned_site_rate = check_pickup_parameter(resDict_SR, "ShortSequence_threshold")
    if re.search("[^\d\.]", Aligned_site_rate):
        print("Error. >ShortSequence_threshold should be floating less than 1.")
        print("Your >Num_rootSequences is", Aligned_site_rate)
        exit()

    dataset = check_pickup_parameter(resDict_SR, "Dataset")
    if dataset == "Exclude3rd" or dataset == "Include3rd" or dataset == "AminoAcid":
        pass
    else:
        print("Error. >Dataset should be Exclude3rd, Include3rd or AminoAcid.")
        print("Your >Dataset is", dataset)
        exit()

    BSthreshold = check_pickup_parameter(resDict_SR, "BSthreshold")
    if re.search("[^\d]", BSthreshold):
        print("Error. >BSthreshold should be integer.")
        print("Your >BSthreshold is", BSthreshold)
        exit()
    elif int(BSthreshold) > 100 or int(BSthreshold) < 50:
        print("Error. >BSthreshold should be less than 100 and more than 50.")
        print("Your >BSthreshold is", BSthreshold)
        exit()

    #treeSearchMethod = check_pickup_parameter(resDict_SR, "TreeSearchMethod")
    treeSearchMethod = "NJ"
    num_rootSequences = check_pickup_parameter(resDict_SR, "Num_rootSequences")
    if re.search("[^\d]", num_rootSequences):
        print("Error. >Num_rootSequences should be integer. 1 is favorable to estimate sister groups.")
        print("Your >Num_rootSequences is", num_rootSequences)
        exit()

    keyNode = check_pickup_parameter(resDict_SR, "KeyNode")
    #name_querySpeciesNode = check_pickup_parameter(resDict_SR, "QuerySpeciesGroup")
    name_querySpecies = check_pickup_parameter(resDict_SR, "QuerySpecies")
    speciesWithGeneFunction = check_pickup_parameter(resDict_SR, "SpeciesWithGeneFunction")
    outdir = check_pickup_parameter(resDict_SR, "Outdir")
    alignment_orthogroups = check_pickup_parameter(resDict_SR, "Alignment_orthogroups")
    mode = check_pickup_parameter(resDict_SR, "Mode")
    if mode == "E" or mode == "D" or mode == "S":
        pass
    else:
        print("Error. Check your >Mode.")
        print(">Mode should be E, D, or S.")
        print("Your >Mode is", mode)
        exit()
    Switch_deleteIntermediateFiles = check_pickup_parameter(resDict_SR, "Switch_deleteIntermediateFiles")
    if Switch_deleteIntermediateFiles == "L" or Switch_deleteIntermediateFiles == "D":
        pass
    else:
        print("Error. Check your mode.")
        print(">Switch_deleteIntermediateFiles should be L or D.")
        print("Your >Switch_deleteIntermediateFiles is", mode)
        exit()
    dbAddress = check_pickup_parameter(resDict_SR, "Database")
    if not os.path.exists(dbAddress):
        print("Error. Your >Database is not found:")
        print(dbAddress)
        exit()
    dbAddress = dbAddress + "/"
    dbLinesTMP, taxonSamplingListTMP = check_pickup_taxonSampling(dbAddress, resDict_SR[">TaxonSampling"])
    #exit("check_pickup_taxonSampling END")
    queryDatabase = ""

    for dbLine in dbLinesTMP:
        if re.search(name_querySpecies + "_", dbLine[0]):
#            queryDatabase = dbLine[2]
            queryDatabase = dbLine[1]
    if not queryDatabase:
        if mode == "E" or mode == "D":
            print("Error: >QuerySpecies", name_querySpecies," cannot be found in >TaxonSampling. Stopped")
            exit()

    return dbLinesTMP, taxonSamplingListTMP, SpeciesTree, blastEvalue, \
Number_of_hits_to_report_per_genome, Aligned_site_rate, dataset, BSthreshold, \
treeSearchMethod, num_rootSequences, keyNode, name_querySpecies, queryDatabase, \
dbAddress, outdir, alignment_orthogroups, mode, Switch_deleteIntermediateFiles, speciesWithGeneFunction


def reorder_dbLines(dbLinesTMP, taxonSamplingListTMP, name_querySpecies, treeFileName):
    #print("#### reorder_dbLines ####")
    #print("name_querySpecies", name_querySpecies)
    leaves = collect_leaves_InOrderFrom_bothNHXnewick(treeFileName)

    # dbLines
    dbLines_treeOrder = []
    for leaf in leaves:
        #print("leaf", leaf)
        for dbLine in dbLinesTMP:
            #print("dbLine[0]", dbLine[0])
            if dbLine[0] == leaf + "_":
               dbLines_treeOrder.append(dbLine)
        #exit()
    #for dbLine in dbLines_treeOrder:
    #    print("dbLine[0]1", dbLine[0])
    #exit()
       
    dbLines = []
    for line in dbLines_treeOrder:
        if not line[0] == name_querySpecies + "_":
            dbLines.append(line)
        #dbLinesTMP1.append(line)
    for line in dbLinesTMP:
        if line[0] == name_querySpecies + "_":
            dbLines.append(line)
    #for dbLine in dbLines:
    #    print("dbLine[0]2", dbLine[0])
    #exit()


    # taxonSamplingList
    taxonSamplingList = []
    for dbline in dbLines:
        for line_ts in taxonSamplingListTMP:
            if re.search(dbline[0], line_ts):
                taxonSamplingList.append(line_ts)
    #for line2 in taxonSamplingList:
    #    print("line2", line2)
    #exit()

    return dbLines, taxonSamplingList


def dirFileMake(queryDatabase, name_querySpecies, queryID):

    recs_queryDB = readFasta_dict(dbAddress, queryDatabase)
    nameline_query = ""
    sequence_query = ""
    for nameline_DB, sequence_DB in recs_queryDB.items():
        #print(nameline_DB)
        if re.search(">" + queryID + " ", nameline_DB) or re.search(">" + queryID + "$", nameline_DB):
            nameline_query = nameline_DB
            sequence_query = sequence_DB
            break
    if not nameline_query:
        print("Error in your control file.")
        print(queryID, "is not found in", dbAddress + queryDatabase)
        exit()

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if not os.path.exists(alignment_orthogroups):
        os.mkdir(alignment_orthogroups)

    if not os.path.exists(eachDirAddress):
        os.mkdir(eachDirAddress)

    fs = open(eachDirAddress + "000_cds_assigned_by_ID.txt", "w")
    nameline_query = re.sub(" .*$", "", nameline_query)
    nameline_query = re.sub(">", ">" + name_querySpecies + "_", nameline_query)
    fs.write(nameline_query + "\n")
    fs.write(sequence_query + "\n")
    fs.close()


def make_treeFile(name_000_uploadedTree, SpeciesTree):
    fhtree = open(eachDirAddress + name_000_uploadedTree, "w")
    fhtree.write(SpeciesTree)
    fhtree.write("\n")
    fhtree.close()


def makeblastdb_database():
    #names_dbFile = files = os.listdir(dbAddress)
    #print("names_dbFile", names_dbFile)
    #exit()
    extensions = ["phr", "pin", "psq"]
    for dbLine in dbLines:
        flag = 0
        for extension in extensions:
            if os.path.isfile(dbAddress + dbLine[1] + "." + extension):
                flag += 1
        if flag < 3:
            print("#### Running makeblastdb.")
            line_makeblastdb = "tools/makeblastdb  -dbtype prot -in  " + dbAddress + dbLine[1]
            #print("line_makeblastdb", line_makeblastdb)
            subprocess.call(line_makeblastdb, shell=True)
### File checking End
##############################################


##############################################
### cDNA and AA file make Start
def checkUplodedFileAsFastaForamt():
    f = open(eachDirAddress + "000_cds_assigned_by_ID.txt")
    lines = list(f)
    f.close()
    if not lines:
        print('Error: Cannot find your sequence file in the specified directory.\n')
        exit()
    if not lines[0].startswith(">"):
        print('Error: Check your sequence file. In fasta format, name line starts with ">"\n')
        exit()
    recs_uploaded = readFasta_dict(eachDirAddress , "000_cds_assigned_by_ID.txt")
#    for name, seq in recs_uploaded.items():
#        if re.search("[^ATGCNX ]", seq):
#            print ('Error: Check your sequence. Sequences should be consist of A,T,G,C,N,X.', "<br>")
#            print ("in ", name, "<br>")
#            exit()


def ckeck_cDNAsequence(recsFN):
    for name, seq in recsFN.items():
        if re.search("[^ATGCNXatgcnx ]", seq):
            print ("Error: Please check sequence in: ", name, "<br>")
            print ("The sequence should be A,T,G,C,N,X. Analysis stopped.")
            exit()


def readPhy_dict(phyFileName):
    phyFile = open(eachDirAddress + phyFileName, "r")
    lines = list(phyFile)
    seqDictFN = OrderedDict()
    for line in lines[1:]:
        name,seq = re.split(" +", line)
        seq = seq.rstrip("\n")
        seqDictFN[">" + name] = seq
    phyFile.close()
    return seqDictFN


def readRes_dict(InfileNameFN):
    #print("InfileNameFN", InfileNameFN)
    flag = 0
    Infile = open(InfileNameFN, "r")
    seqDictFN  = OrderedDict()
    stock = []
    for Line in Infile:
        #print("Line", Line)
        if re.search("^#", Line) or re.search("^$", Line):
            continue
        Line = Line.rstrip("\n")
        if Line[0] == ">":
            if flag == 0:
                Name = Line
                flag = 1
            else:
                Name = re.sub(" +$", "", Name)
                seqDictFN[Name] = stock
                Name = Line
                stock = []
        else:
            stock.append(Line)
    Name = re.sub(" +$", "", Name)
    seqDictFN[Name] = stock
    Infile.close()

    return seqDictFN


def readFasta_dict(dirAddressFN, InfileNameFN):
    #print("dirAddressFN + InfileNameFN", dirAddressFN + InfileNameFN)
    #exit()
    Infile = open(dirAddressFN + InfileNameFN, "r")
    seqDictFN  = OrderedDict()
    for Line in Infile:
        Line = Line.rstrip("\n")
        if not Line:
            continue
        elif Line[0] == ">":
            Name            = Line
            seqDictFN[Name] = ""
        else:
            Line = Line.replace("\n", "")
            Line = Line.replace("\r", "")
            #if InfileNameFN == "000_cds_assigned_by_ID.txt":
            #    Line = re.sub("-", "", Line)
            seqDictFN[Name] += Line.upper()
    Infile.close()
    return seqDictFN


def splitDna(dna):
    codons = []
    for start in range(0, len(dna)-2, 3):
        codons.append(dna[start:start+3])
    return(codons)


def translation(dna):
    dna = dna.upper()
    protein = ""
    for codon in splitDna(dna):
        aa = geneticCode.get(codon, "X")
        protein = protein + aa
    return protein


def aaSeqMaker():
    recfn = readFasta_dict(eachDirAddress , "000_cds_assigned_by_ID.txt")

    fa = open(eachDirAddress + "000_translated_cds_assigned_by_ID.txt","w")

    # >Human_ENSP00000259365 gene:ENSG00000136842 transcript:ENST00000259365 gene_biotype:protein_coding
    for name, seq in recfn.items():
        newName = re.sub("\|", "", name)
        fa.write(newName + "\n")
#        fa.write(translation(seq) + "\n")
        fa.write(seq + "\n")
    fa.close()


### cDNA and AA file make End
##############################################


###############################################################
##################### Tree manipulation START
def collect_leaves_InOrderFrom_bothNHXnewick(treeFileName):
    #print("### collect_leaves_InOrderFrom_bothNHXnewick ###")
    #print("treeFileName", treeFileName)
    #exit()
    treeTMP = open(eachDirAddress + treeFileName, "r")
    tree = list(treeTMP)[0]
    treeTMP.close()
    leaves = tree.split(",")
    leaves = [re.sub("^\(*([^:\(\)]+).*$", r"\1", leaf) for leaf in leaves]
    leaves = [re.sub("[\t\n]", "", leaf) for leaf in leaves]
    #for leaf in leaves:
    #    print(leaf, "<br>")
    #exit()
    return leaves


def test_querySpecies_is_in_speciesTree(taxonSamplingList, SpeciesTree):
    speciesNames_taxonSamplingList = [] 
    for speciesName_taxonSamplingList in taxonSamplingList:
        speciesName_taxonSamplingList = re.sub("_[^_]+$", "", speciesName_taxonSamplingList)
        speciesNames_taxonSamplingList.append(speciesName_taxonSamplingList)

    if re.search("\)\)", SpeciesTree) or re.search("\),", SpeciesTree) :
        print ("Error: In the species tree, all nodes should have node name.")
        print ("Please check you newick format using FigTree or TreeGraph_2.")
        exit()

    num_rightOpen = SpeciesTree.count("(")
    num_leftOpen  = SpeciesTree.count(")")
    num_comma     = SpeciesTree.count(",")
    if num_rightOpen != num_leftOpen:
        print ("Error: In the species tree, numbers of right and left parehtheses should be same. <br>")
        exit()
    if num_rightOpen != num_comma or num_leftOpen != num_comma :
        print ("Error: In the species tree, all nodes should be bufuricated.")
        exit()

    speciesNames_in_speciesTree = ""
    for node in allNodes_speciesTree:
        #print("node", node)
        if re.search("S=" + keyNode + ":", node[2]):
            speciesNames_in_speciesTree = node[1]
    if not speciesNames_in_speciesTree:
        print("Error: KeyNode, ", keyNode, ", is not found in the SpeciesTree:")
        exit()

    allSpeciesNames_speciesTree = allNodes_speciesTree[0][1]
    for speciesName_uploadedSeqFile in speciesNames_taxonSamplingList:
        #print("speciesName_uploadedSeqFile", speciesName_uploadedSeqFile)
        if not speciesName_uploadedSeqFile in allSpeciesNames_speciesTree:
            print("Error: Sequence, ", speciesName_uploadedSeqFile, ", is not found in your species tree:")
            print(SpeciesTree)
            exit()
    for dbLine in dbLines:
        spName_dbLine = dbLine[0]
        spName_dbLine = re.sub("_$","",spName_dbLine)
        if not spName_dbLine in allSpeciesNames_speciesTree:
            print ("Error: ", spName_dbLine, ", is not found in your uploaded-species tree.")
            exit()


def change_nhx_to_newick_with_NHXnodeName(tree_NHX):
  
    cladeReg = "\)([^\[]+\[&&NHX[^\]]+\])"
    count = 0
    while re.search(cladeReg, tree_NHX):
        #print("count   :", count)
        matchA = re.search(cladeReg, tree_NHX)
        expTemp  = matchA.group(1);
        #print("expTemp:", expTemp)
        #if count == 10:
        #    exit()

        matchB = re.search("(\[.*\])", expTemp)
        exp  = matchB.group(1)
 
        if re.search("B=([\d]+)", expTemp):
            matchC = re.search("B=([\d]+)", expTemp)
            bs = matchC.group(1)
        else:
            bs = "r"

        exp = re.sub("&&NHX:", "", exp)
        exp = re.sub(":",      "_", exp)
        exp = re.sub("_B=.*$", "", exp)
        
        #print("expTemp :", expTemp)
        #print("bs      :", bs)
        #print("exp     :", exp)
        #print("bs + exp:", bs + exp)

        tree_NHX = re.sub(cladeReg, ")" + bs + "_" + exp, tree_NHX, 1)
        count += 1
        #print("tree2: ", tree_NHX)
        #print()
    
    tree_NHX = re.sub("\:[\d|\.E-]+", "", tree_NHX)   #Delete the branch lengths

    ## left node name
    tree_NHX = re.sub("\[&&NHX:[^]]+\]", "",  tree_NHX)
    tree_NHX = re.sub("[\[\]]",          "", tree_NHX)
    #tree_NHX = re.sub("\]",              "",  tree_NHX)

    return tree_NHX;


def collect_nodes_from_speciesTree():
    tree_newick = SpeciesTree
    nodes = []     # 2D array for nodes
    cladeReg = "\(([^\(\)]+)\)(\w+)"
    while re.search(cladeReg, tree_newick):
        tree_newick = tree_newick.rstrip("\n")
        match = re.search(cladeReg, tree_newick)             # Pick up the smallest and leftmost clade for the following analysis
        leavesString = match.group(1)
        exp = match.group(2)
        #leaves = set(leavesString.split(","))
        leaves = leavesString.split(",")
        tree_newick = re.sub(cladeReg, r"\1", tree_newick, 1)      # Delete the outer parentheses from the analyzed clade
        nodes.append([len(leaves), leaves, "[&&NHX:S=" + exp + ":D=N:B=speciesTree]"])

    sortedNodes     = sorted(nodes, key=lambda x:x[0], reverse=True)

    ## Add leaf as a clade
    
    largestClade = sortedNodes[0];
    for leaf in largestClade[1]:
        exp           = leaf
        #tempLeafClade = set([leaf])
        sortedNodes.append([1, [leaf], "[&&NHX:S=" + exp + "]"])

    #for node in sortedNodes:
    #    print(node[0])
    #    print(node[1])
    #    print(node[2])
    #    print("")
    #exit()
    return sortedNodes


'''
def collect_nodes_from_newick(tree_newick):

    nodes = []     # 2D array for nodes
    cladeReg = "\(([^\(\)]+)\)(\w+)"
    while re.search(cladeReg, tree_newick):
        tree_newick = tree_newick.rstrip("\n")
        match = re.search(cladeReg, tree_newick)             # Pick up the smallest and leftmost clade for the following analysis
        leavesString = match.group(1)
        exp = match.group(2)
        #leaves = set(leavesString.split(","))
        leaves = leavesString.split(",")
        #print("leaves", leaves)
        #exit()
        tree_newick = re.sub(cladeReg, r"\1", tree_newick, 1)      # Delete the outer parentheses from the analyzed clade
        nodes.append([len(leaves), leaves, "[&&NHX:S=" + exp + ":D=N:B=speciesTree]"])

    sortedNodes     = sorted(nodes, key=lambda x:x[0], reverse=True)

    ## Add leaf as a clade
    
    largestClade = sortedNodes[0];
    for leaf in largestClade[1]:
        exp           = leaf
        tempLeafClade = set([leaf])
        sortedNodes.append([1, tempLeafClade, "[&&NHX:S=" + exp + "]"])

    return sortedNodes
'''


def collect_nodes_from_NHX(treeFN):
    clades = []     # 2D array for clades
    expReg = "\[.*?\]"
    cladeReg = "\(([^\(\)]+)\)(.*?\[.*?\])"
    while re.search(cladeReg, treeFN):
        treeFN = treeFN.rstrip("\n")
        match = re.search(cladeReg, treeFN)                # Pick up the smallest and leftmost clade for the following analysis
        leavesString = match.group(1)
        exp = match.group(2)
        leavesString = re.sub("\[.*?\]", "", leavesString)   # Delete exp [...] of internal branches
        leavesString = re.sub(":\d+\.\d+E-\d*", "", leavesString)   # Delete blanch lengths with E-
        leavesString = re.sub(":\d+\.\d+", "", leavesString)       # Delete blanch lengths of internal branches
        #leaves = set(leavesString.split(","))
        leaves = leavesString.split(",")
        treeFN = re.sub(cladeReg, r"\1", treeFN, 1)                # Delete the outer parentheses from the analyzed clade
        clades.append([len(leaves), leaves, exp])

    sortedClades = sorted(clades, key=lambda x:x[0], reverse=True)

    ## Add leaf as a clade
    largestClade = sortedClades[0];
    for leaf in largestClade[1]:
        #print ("leaf:", leaf)
        match         = re.search(r"([^_]+)_", leaf)
        exp           = "[&&NHX:S=" + match.group(1) +"]"
        #tempLeafClade = set([leaf])
        #tempLeafClade = set([leaf])
        sortedClades.append([1, [leaf], exp])

    return sortedClades


def identify_orthogroup(treeNHX):
    #print("### identify_orthogroup ###")
    #print("treeNHX", treeNHX)
    #f = open(eachDirAddress + "000_speciesTree.txt")
    #speciesTree = "".join(list(f))
    #speciesTree = re.sub("[ \n]", "", speciesTree)
    #f.close()

    allGeneNodesSR = collect_nodes_from_NHX(treeNHX)

    topHits = topHitPicker()
    #print("topHits", topHits)
    #exit()
    nameLine_cds_blastTopHit = list(topHits.values())[0][0]
    nameLine_cds_assigned_by_ID = list(topHits.keys())[0]
    
    #print("nameLine_cds_blastTopHit", nameLine_cds_blastTopHit)
    #exit()
    #print("nameLine_cds_assigned_by_ID", nameLine_cds_assigned_by_ID)
    #print("keyNode", keyNode)
    
    candidates_orthogroup = []
    #print("nameLine_cds_blastTopHit[1:]", nameLine_cds_blastTopHit[1:])
    for node in allGeneNodesSR:
        #print("node", node)
        #exit()        
        criterion = 0
        for leaf in node[1]:
            #print("leaf", leaf)
            if leaf == nameLine_cds_blastTopHit[1:]:
            #if re.search(nameLine_cds_blastTopHit[1:], leaf):
                criterion += 1


        if re.search("S=" + keyNode, node[2]):
            criterion += 1

        if criterion == 2:
           candidates_orthogroup.append(node)

    #print("")

    orthogroupSR = []
    if not candidates_orthogroup:
        
        flag = 0
        leaves = leafCollectInOrderFrom_bothNHXnewick(treeNHX)
        for leaf in leaves:
            if leaf == nameLine_cds_blastTopHit[1:]:
                flag = 1
                break
        if flag == 1:
            orthogroupSR = [0, 0, "noOrthogroup_noKeynode"]
        else:
            orthogroupSR = [0, 0, "noOrthogroup_noQuerySequence"]
    else:
        orthogroupSR = candidates_orthogroup.pop()
        #print("orthogroupSR", orthogroupSR)

    return orthogroupSR


def leafCollectInOrderFrom_bothNHXnewick (tree):
    leaves = tree.split(",")
    leaves = [re.sub("^\(*([^:\(\)]+).*$", r"\1", leaf) for leaf in leaves]
    leaves = [re.sub("[\t\n]", "", leaf) for leaf in leaves]
    #for leaf in leaves:
    #    print(leaf, "<br>")
    #exit()
    return leaves


def collect_speciesNames_in_orthogroup():
    #allNodes_speciesTree = collect_nodes_from_newick(SpeciesTree)
    leaves_species = leafCollectInOrderFrom_bothNHXnewick(SpeciesTree)
    orthoSpeciesGroup = ""
    speciesNames_in_orthogroup_FN = []
    for node in allNodes_speciesTree:
        if re.search("S=" + keyNode, node[2]):
            orthoSpeciesGroup = node
            break
    for leaf_species in leaves_species:
        if leaf_species in orthoSpeciesGroup[1]:
            speciesNames_in_orthogroup_FN.append(leaf_species)
    return speciesNames_in_orthogroup_FN


def identify_speciesNode(name_node):
    for node in allNodes_speciesTree:
        if re.search("S=" + name_node, node[2]):
            return node


def identify_targetGeneNode(allNodes_SR, name_speciesNode, separationType, queryGeneLeaf):
    #print("### Start identify_targetGeneNode")
    #print("name_speciesNode", name_speciesNode)
    #print("separationType", separationType)
    #targetSpeciesNode = identify_speciesNode(name_speciesNode)
    candidates_queryGeneNode = []

    flag_counting = 0
    for node in allNodes_SR:
        flag_counting += 1
        if flag_counting == 1:
            continue

        if separationType == "SisterGeneGroups" and flag_counting == 2:
            continue

        #print("node[2]", node[2])


        #print("node[2]", node[2])

        criterion = 0

        for leaf in node[1]:
            if leaf == queryGeneLeaf:
                criterion += 1

        keyWord = "S=" + name_speciesNode + "[:\]]"
        if re.search(keyWord, node[2]):
            #print("keyWord2", keyWord)
            criterion += 1

        if criterion == 2:
           candidates_queryGeneNode.append(node)
           
        #print("")

    node_identified = []
    if not candidates_queryGeneNode:
        node_identified = [0, 0, "NoGeneNode"]
    else:
        #node_identified = candidates_queryGeneNode.pop()
        node_identified = candidates_queryGeneNode[0]

    #print("node_identified", node_identified)
    #print("### END\n")
    return node_identified


#def identify_focalGeneNode_with_focalNodeName_and_parentNodeName(allNodesSR, names_speciesNodes, name_queryLeaf):
#    keyNodes   = []
#    for clade in allNodesSR:
#        criterion_focalSpeciesClade = 0
#        
#        for name_focalGeneNode in names_speciesNodes:
#            if re.search("S=" + name_focalGeneNode + ":", clade[2]):
#                #print("Hit1")
#                criterion_focalSpeciesClade += 1
#
#        if name_queryLeaf[1:] in clade[1]:
#            #print("Hit2")
#            criterion_focalSpeciesClade += 1
#
#        if criterion_focalSpeciesClade == 2:
#            keyNodes.append(clade)
#
#    if not keyNodes:
#        print ("Error: Stopped in species tree (Newick). No focal gene clade was found for ", names_speciesNodes, querySpeciesName)
#        exit()
#    else:
#        focalNode = keyNodes.pop()
#        return focalNode


def isIndependent(checkCladeLeaves, focalCladeLeaves):
    for ckeckCladeLeaf in checkCladeLeaves:
        #if [focalCladeLeaf for focalCladeLeaf in focalCladeLeaves if re.search(ckeckCladeLeaf, focalCladeLeaf)]:
        for focalCladeLeaf in focalCladeLeaves:
            if ckeckCladeLeaf == focalCladeLeaf:
                return 0
    return 1;


def isOverlapped(checkCladeLeaves, focalCladeLeaves):
    for ckeckCladeLeaf in checkCladeLeaves:
        #if [focalCladeLeaf for focalCladeLeaf in focalCladeLeaves if re.search(ckeckCladeLeaf, focalCladeLeaf)]:
        for focalCladeLeaf in focalCladeLeaves:
            if ckeckCladeLeaf == focalCladeLeaf:
                return 1
    return 0;


def count_species_in_gene_clade(checkSpeciesCladeLeaves, focalGeneCladeLeaves):
    hits_sr = 0
    for ckeckSpeciesCladeLeaf in checkSpeciesCladeLeaves:
        hits_4_each_species = 0
        for focalGeneCladeLeaf in focalGeneCladeLeaves:
            if re.search(ckeckSpeciesCladeLeaf, focalGeneCladeLeaf):
                hits_4_each_species += 1
        if hits_4_each_species > 0: 
            hits_sr += 1
    return hits_sr


def collect_childNodesincluding_querySpecies(allGeneNodesSR, focalNode_SR):
    childSpeciesNodes_orthogorup = collect_childNodes(allGeneNodesSR, focalNode_SR)
    childNodesincluding_querySpecies = []
    for speciesNode in childSpeciesNodes_orthogorup:
        if name_querySpecies in speciesNode[1]:
            childNodesincluding_querySpecies.append(speciesNode)
    return childNodesincluding_querySpecies


def collect_childNodes(allGeneNodesSR, focalNode_SR):
    #print("collect_childNodes, focalNode_SR", focalNode_SR)
    childNodesSR = []
    for each_Node in allGeneNodesSR:
        leaves_each_Node = set(each_Node[1])
        leaves_focalNode_SR = set(focalNode_SR[1])
        #print("leaves_each_Node", leaves_each_Node)
        #print("leaves_focalNode_SR", leaves_focalNode_SR)
        #if eachNode[1].issubset(focalNodeSR[1]):
        if leaves_each_Node.issubset(leaves_focalNode_SR):
            childNodesSR.append(each_Node)
        #print("")
    #print("len(childNodesSR)", len(childNodesSR))
    #for node in childNodesSR:
    #    print("node[2]", node[2])
    #print("")
    return childNodesSR


def collect_childBranchLavels(SpeciesTree, nodeName, name_querySpecies):
    #print("collect_childBranchLavels")
    #print("nodeName", nodeName)
    #print("name_querySpecies", name_querySpecies)
    #allNodes_speciesTree  = collect_nodes_from_newick(SpeciesTree)
    flag = 0

    hildBranchLavels = []
    for node_speciesTree in allNodes_speciesTree:
        branchLabel = node_speciesTree[2]
        #print("branchLabel", branchLabel)
        if flag == 1:
            if name_querySpecies in node_speciesTree[1]:
                match = re.search("S=([^:\]]+)[:\]]", branchLabel)
                hildBranchLavels.append(match.group(1))
        if re.search("S=" + nodeName + ":", node_speciesTree[2]):
            hildBranchLavels.append(nodeName)
            flag = 1

    return hildBranchLavels


def collect_sisterGroups(allGeneNode_SR, focalNode_SR):
    #print("## collect_sisterGroups ##")
    #print("focalNode_SR", focalNode_SR)
    sisterGeneGroups_gettingDeeper_SR = []
    ancestralNodesDecrement2 = collect_ancestralNodes(allGeneNode_SR, focalNode_SR)
    #print("ancestralNodesDecrement2")
    #for node in ancestralNodesDecrement2:
    #    print("node", node)
    #exit()
    for ancestralNode in reversed(ancestralNodesDecrement2):
        if ancestralNode == focalNode_SR:
            #print("ssddfdd")
            continue
        daughterNode_1st, daughterNode_2nd = identify_daughterNodes(allGeneNode_SR, ancestralNode)
        leaves_daughterNode_1st = set(daughterNode_1st[1])
        leaves_focalNode_SR = set(focalNode_SR[1])
        #if isOverlapped(daughterNode_1st[1], focalNode_SR[1]):
        if leaves_focalNode_SR.issubset(leaves_daughterNode_1st):
            sisterGeneGroups_gettingDeeper_SR.append(daughterNode_2nd)
        else:
            sisterGeneGroups_gettingDeeper_SR.append(daughterNode_1st)
    #print("sisterGeneGroups_gettingDeeper_SR", sisterGeneGroups_gettingDeeper_SR)
    #exit()
    return sisterGeneGroups_gettingDeeper_SR


def collect_ancestralNodes(allNodes_SR, node_SR):
    ancestralGroupsDecrementSR = []
    for eachLargerNode in allNodes_SR:
        leaves_node_SR = set(node_SR[1])
        leaves_eachLargerNode = set(eachLargerNode[1])
        if leaves_node_SR.issubset(leaves_eachLargerNode):
            ancestralGroupsDecrementSR.append(eachLargerNode)
    return ancestralGroupsDecrementSR


def identify_parentNode(ancestralNodesDecrementSR, ancestralDepthSR, focalNodeSR):
    parentalNodeSR = ""
    if ancestralNodesDecrementSR[0][2] == focalNodeSR[2]:
        return [0, 0, "NoParentalGeneNode"]
    ancestralNodesSRIncrement = ancestralNodesDecrementSR[::-1]
    for i in range(ancestralDepthSR, -1, -1):
        if ancestralNodesSRIncrement[i]:
            parentalNodeSR = ancestralNodesSRIncrement[i]
            return parentalNodeSR


#def identify_ancestralSpeciesNode(filename_uploadedSpciesTree):
#    focalSpeciesNode = identify_node(filename_uploadedSpciesTree)
#    #allNodes_speciesTree = collect_nodes_from_newick(SpeciesTree)
#    ancestralDepth = 1
#    ancestralSpeciesNodesDecrement = collect_ancestralNodes(allNodes_speciesTree, focalSpeciesNode)
#    parentNode = identify_parentNode(ancestralSpeciesNodesDecrement, ancestralDepth, keyNode)
#    match = re.search("S=([^:]+):", parentNode[2])
#    name_ancestralSpeciesNode = match.group(1)
#    return name_ancestralSpeciesNode


def identify_sisterNode(childNodes_parentNode_SR, focalNode_SR):
    for eachChildNode in childNodes_parentNode_SR:
        if isIndependent(eachChildNode[1], focalNode_SR[1]):
            return eachChildNode;


def identify_daughterNodes(allNodesSR, focalNodeSR):
    childNodes = collect_childNodes(allNodesSR, focalNodeSR)
    daughterNode1st = []
    daughterNode2nd = []
    if childNodes:
        daughterNode1st = childNodes[1];
        daughterNode2nd = identify_sisterNode(childNodes, daughterNode1st)
    return(daughterNode1st, daughterNode2nd)


def identify_sisterGeneNode(allGeneNodes_SR, targetGeneNode_FN):
    parentNode_of_targetGeneNode_SR = ""
    sisterGeneNode_SR = ""
    ancestralDepth = 1
    ancestralNodesDecrement = collect_ancestralNodes(allGeneNodes_SR, targetGeneNode_FN)
    parentNode_of_targetGeneNode_SR = identify_parentNode(ancestralNodesDecrement, ancestralDepth, targetGeneNode_FN)
    if parentNode_of_targetGeneNode_SR[2] == "NoParentalGeneNode":
        sisterGeneNode_SR = [0, 0, "NoSisterNode"]
    else:
        childNodes_parentGeneNode = collect_childNodes(allGeneNodes_SR, parentNode_of_targetGeneNode_SR)
        sisterGeneNode_SR = identify_sisterNode(childNodes_parentGeneNode, targetGeneNode_FN)
    return parentNode_of_targetGeneNode_SR, sisterGeneNode_SR


def identifiy_orthogroup_speciesNode():
    orthogroup_speciesNode = ""
    for node in allNodes_speciesTree:
        if re.search("S=" + keyNode + ":", node[2]):
            orthogroup_speciesNode = node
    return orthogroup_speciesNode


def count_duplications_for_speciesNodes(allGeneNodesSR, topHitName_1stQuery):
    #print("topHitName_1stQuery", topHitName_1stQuery)
    targetGeneNode = ""
    for node in allGeneNodesSR:
        if node[0] == 1 and node[1][0] == topHitName_1stQuery:
            targetGeneNode = node

    geneNodes_including_querySequence = collect_ancestralNodes(allGeneNodesSR, targetGeneNode)
    #for node in geneNodes_including_querySequence:
    #    print(node[2])
    #exit()

    recs_duplications_for_speciesNodes_FN = OrderedDict()
    #for childSpeciesNode_of_orthogroup in childSpeciesNodes_orthogorup:
    #    if not name_querySpecies in childSpeciesNode_of_orthogroup[1]:
    #        continue
    for childSpeciesNode_of_orthogroup in childSpeciesNodes_orthogroup_including_querySpecies:

        speciesNodeName = make_nodeName_from_nodeLavel_NHXstyle(childSpeciesNode_of_orthogroup[2])
        #print("speciesNodeName", speciesNodeName)
        count_dup = 0
        flag = 0
        for node_geneTree in geneNodes_including_querySequence:
            #print("node_geneTree[2]", node_geneTree[0], node_geneTree[2])
            if flag == 1:
                if re.search(":D=Y[:\]]", node_geneTree[2]):
                    match = re.search("^([^\[]+)\[.*S=([^:]+):", node_geneTree[2])
                    nodeID_geneTree = match.group(1)
                    geneNodeName = match.group(2)
                    #print("nodeID_geneTree", nodeID_geneTree)
                    #print("geneNodeName", geneNodeName)
                    if geneNodeName == speciesNodeName and nodeID_geneTree.startswith("n"):
                        count_dup += 1
                        #print("count")

            #if flag == 1:
            #    if re.search("S=" + speciesNodeName + ":D=Y", node_geneTree[2]):
            #        count_dup += 1

            flag = 1
        #print("")

        recs_duplications_for_speciesNodes_FN[speciesNodeName] = count_dup

    return recs_duplications_for_speciesNodes_FN
########################################################################################################################################


##################### Tree manipulation END
###############################################################


###############################################################
### makeSummary Start
def count_blastHits(infile):
    f = open(eachDirAddress + infile)
    lines = list(f)
    f.close()
    
    spNamePrefixes = []
    for spNameTMP in taxonSamplingList:
        match = re.search("^([^_]+_)([^_]+)$", spNameTMP)
        spName = match.group(1)
        spNamePrefixes.append(spName)

    rec_blastHitsNums = OrderedDict()
    for spNamePrefix in spNamePrefixes:
        hits = len([line for line in lines if re.search(spNamePrefix, line)])
        rec_blastHitsNums[spNamePrefix[:-1]]= hits
    return rec_blastHitsNums


def error_makeSummary(resultFN):
    fs = open(eachDirAddress + "100_analysisSummary.txt", "w")

    fs.write("################ Results: 1st analysis ################\n\n")
    fs.write(">BS_of_orthogroupBasalNode\n")
    fs.write(resultFN + "\n");
    fs.write("\n")

    recs_cds_assigned_by_ID = readFasta_dict(eachDirAddress, "000_cds_assigned_by_ID.txt")
    fs.write(">QuerySequence\n")
    lines_hit_query = make_lines_hit_query(recs_cds_assigned_by_ID)
    for line in lines_hit_query:
        fs.write(line)
    fs.write("\n")

    fs.close()


def make_bsvalue_from_nodeLavel_NHXstyle(nodeLavel_NHXstyle):
    #print("nodeLavel_NHXstyle", nodeLavel_NHXstyle)
    if nodeLavel_NHXstyle.startswith("["):
        return("leaf")
    match = re.search("B=([\d]+)", nodeLavel_NHXstyle)
    if match:
        return(match.group(1))
    else:
        return("r")

def make_duplicationStatus_from_nodeLavel_NHXstyle(nodeLavel_NHXstyle):
    #print("nodeLavel_NHXstyle", nodeLavel_NHXstyle)
    if nodeLavel_NHXstyle.startswith("["):
        return("leaf")
    if re.search(":D=Y[:\]]", nodeLavel_NHXstyle):
        return("D=Y")
    if re.search(":D=N[:\]]", nodeLavel_NHXstyle):
        return("D=N")
    else:
        print("Error in  make_duplicationStatus_from_nodeLavel_NHXstyle.")
        print("Check nodeLavel", nodeLavel_NHXstyle)
        exit()


def make_nodeName_from_nodeLavel_NHXstyle(nodeLavel_NHXstyle):
    match = re.search("S=([^:]+)[:\]]", nodeLavel_NHXstyle)
    if match:
        return(match.group(1))
    else:
        return("No_node_lavel")


def check_sisterGroupName_included_in_ancestralSpeciesNodeNames(allGeneNodesSR, targetSpeciesNode, sisterGeneGroup_SR):
    #targetSpeciesNode = [8, ['Scleropages-formosus', 'Paramormyrops-kingsleyae', 'Megalops-cyprinoides', 'Anguilla-anguilla', 'Clupea-harengus', 'Danio-rerio', 'Gasterosteus-aculeatus', 'Oryzias-latipes'], '[&&NHX:S=Teleostei:D=N:B=speciesTree]']
    #sisterGeneGroup_SR = [3, ['Erpetoichthys-calabaricus_XP028657515.1_protocadherin-18-is', 'Acipenser-ruthenus_XP033890080.2_protocadherin-18a-isoform-', 'Lepisosteus-oculatus_ENSLOCT00000013009.1_pcdh18b-protocadhe'], 'n33:0.008780038052[&&NHX:S=Actinopterygii:D=Y:B=74.0]']
    ##print("###########")
    #print("targetSpeciesNode", targetSpeciesNode)
    #print("sisterGeneGroup_SR", sisterGeneGroup_SR)
    #print("")

    nodeName_sisterGeneGroup_SR = make_nodeName_from_nodeLavel_NHXstyle(sisterGeneGroup_SR[2])
    #print("nodeName_sisterGeneGroup_SR", nodeName_sisterGeneGroup_SR)
    if sisterGeneGroup_SR[0] == 1:
        return nodeName_sisterGeneGroup_SR

    names_ancestralGeneNode = []
    ancestralGeneNodes_of_querySpeciesNode = collect_ancestralNodes(allNodes_speciesTree, targetSpeciesNode)
    for ancestralGeneNode_of_querySpeciesNode in ancestralGeneNodes_of_querySpeciesNode:
        names_ancestralGeneNode.append(make_nodeName_from_nodeLavel_NHXstyle(ancestralGeneNode_of_querySpeciesNode[2]))

    if nodeName_sisterGeneGroup_SR in names_ancestralGeneNode:
        names_childGeneNode = []
        '''
        childNodes = collect_childNodes(allGeneNodesSR, sisterGeneGroup_SR)
        for childGeneNode in childNodes:

            if childGeneNode[0] == 1:
                continue

            #print("childGeneNode[2]", childGeneNode[2])
            name_childGeneNode = make_nodeName_from_nodeLavel_NHXstyle(childGeneNode[2])

            ### contain node names
            #if name_childGeneNode != nodeName_sisterGeneGroup_SR:
            if name_childGeneNode != nodeName_sisterGeneGroup_SR and name_childGeneNode not in names_childGeneNode:
                names_childGeneNode.append(name_childGeneNode)
            #print("names_childGeneNode1", names_childGeneNode)

            ### contain one-leaf clade
            #daughterGeneNode_1st, daughterGeneNode_2nd = identify_daughterNodes(allGeneNodesSR, childGeneNode)
            #if daughterGeneNode_1st[0] == 1 and daughterGeneNode_2nd[0] > 1:
            #    names_childGeneNode.append(make_nodeName_from_nodeLavel_NHXstyle(daughterGeneNode_1st[2]))
            #if daughterGeneNode_2nd[0] == 1 and daughterGeneNode_1st[0] > 1:
            #    names_childGeneNode.append(make_nodeName_from_nodeLavel_NHXstyle(daughterGeneNode_2nd[2]))

            print("names_childGeneNode2", names_childGeneNode)
            print("")

        ### contain one-leaf clade derived from basal separation in the sisterGeneGroup_SR
        #daughterGeneNode_1st, daughterGeneNode_2nd = identify_daughterNodes(allGeneNodesSR, sisterGeneGroup_SR)
        #if daughterGeneNode_1st[0] == 1:
        #    names_childGeneNode.append(make_nodeName_from_nodeLavel_NHXstyle(daughterGeneNode_1st[2]))
        #if daughterGeneNode_2nd[0] == 1:
        #    names_childGeneNode.append(make_nodeName_from_nodeLavel_NHXstyle(daughterGeneNode_2nd[2]))
        '''

        ### contain daughter node clades derived from basal separation in the sisterGeneGroup_SR
        daughterGeneNode_1st, daughterGeneNode_2nd = identify_daughterNodes(allGeneNodesSR, sisterGeneGroup_SR)

        candidateName_sisterGeneNode_1stDaugter = make_nodeName_from_nodeLavel_NHXstyle(daughterGeneNode_1st[2])
        #print("daughterGeneNode_1st", daughterGeneNode_1st)
        #print("names_ancestralGeneNode", names_ancestralGeneNode)
        if daughterGeneNode_1st[0] == 1 or candidateName_sisterGeneNode_1stDaugter not in names_ancestralGeneNode:
            names_childGeneNode.append(candidateName_sisterGeneNode_1stDaugter)
        else:
            daughterGeneNode1_1st, daughterGeneNode1_2nd = identify_daughterNodes(allGeneNodesSR, daughterGeneNode_1st)
            temp1 = make_nodeName_from_nodeLavel_NHXstyle(daughterGeneNode1_1st[2])
            temp2 = make_nodeName_from_nodeLavel_NHXstyle(daughterGeneNode1_2nd[2])
            names_childGeneNode.append(temp1)
            names_childGeneNode.append(temp2)

        candidateName_sisterGeneNode_2ndDaugter = make_nodeName_from_nodeLavel_NHXstyle(daughterGeneNode_2nd[2])
        if daughterGeneNode_2nd[0] == 1 or candidateName_sisterGeneNode_2ndDaugter not in names_ancestralGeneNode:
            names_childGeneNode.append(candidateName_sisterGeneNode_2ndDaugter)
        else:
            daughterGeneNode2_1st, daughterGeneNode2_2nd = identify_daughterNodes(allGeneNodesSR, daughterGeneNode_2nd)
            temp1 = make_nodeName_from_nodeLavel_NHXstyle(daughterGeneNode2_1st[2])
            temp2 = make_nodeName_from_nodeLavel_NHXstyle(daughterGeneNode2_2nd[2])
            names_childGeneNode.append(temp1)
            names_childGeneNode.append(temp2)

        unique_names_childGeneNode = sorted(set(names_childGeneNode), key=names_childGeneNode.index)
        nodeName_sisterGeneGroup_SR = nodeName_sisterGeneGroup_SR + "(" + ";".join(unique_names_childGeneNode) + ")"

    #print("nodeName_sisterGeneGroup_SR", nodeName_sisterGeneGroup_SR)
    #print("")
    #exit()
    
    return nodeName_sisterGeneGroup_SR


def add_makeSummary(outfile_summary2):
    #print("#### add_makeSummary ####")
    fSum = open(eachDirAddress + outfile_summary2, "a")

    resDict_1stSummary = readRes_dict(eachDirAddress + "100_analysisSummary.txt")
    topHitName_1stQuery = resDict_1stSummary[">QuerySequence"][0]
    topHitName_1stQuery = re.sub(" +.*", "", topHitName_1stQuery)
    #print("topHitName_1stQuery", topHitName_1stQuery)
    #exit()
    #topHitName_1stQuery = re.sub(" .*$", "", topHitName_1stQuery)

    fSum.write("################ Results: 2nd analysis ################\n\n")

    secondRearrangedTreeTMP = open(eachDirAddress + "230_2ndtreeRootBS100.txt.rearrange.0")
    rearranged_2nd_gene_tree_NHX = list(secondRearrangedTreeTMP)[0]
    secondRearrangedTreeTMP.close()

    fSum.write(">TreeSearchMethod\n")
    fSum.write(treeSearchMethod)
    fSum.write("\n\n")

    fSum.write(">2nd_rearranged_gene_tree_newick\n")
    rearrangedTreeNewick = change_nhx_to_newick_with_NHXnodeName(rearranged_2nd_gene_tree_NHX)
    fSum.write(rearrangedTreeNewick)
    fSum.write("\n")

    fSum.write(">2nd_rearranged_gene_tree_NHX\n")
    fSum.write(rearranged_2nd_gene_tree_NHX)
    fSum.write("\n")

    secondTreeTMP = open(eachDirAddress + "230_2ndtree.txt")
    secondTree = list(secondTreeTMP)[0]
    secondTreeTMP.close()
    fSum.write(">2nd_gene_tree_newick\n")
    fSum.write(secondTree)
    fSum.write("\n")

    allGeneNodesSR = collect_nodes_from_NHX(rearranged_2nd_gene_tree_NHX)

    fSum.write(">MonophyleticGeneGroups\n")
    #for childSpeciesNode_of_orthogroup in childSpeciesNodes_orthogorup:
    #    if not name_querySpecies in childSpeciesNode_of_orthogroup[1]:
    #        continue
    for childSpeciesNode_of_orthogroup in childSpeciesNodes_orthogroup_including_querySpecies:

        name_childSpeciesNode_of_orthogroup = make_nodeName_from_nodeLavel_NHXstyle(childSpeciesNode_of_orthogroup[2])
        #print("name_childSpeciesNode_of_orthogroup1", name_childSpeciesNode_of_orthogroup)
        targetGeneNode = identify_targetGeneNode(allGeneNodesSR, name_childSpeciesNode_of_orthogroup, "MonophyleticGeneGroups", topHitName_1stQuery)
        #print("targetGeneNode", targetGeneNode[0], targetGeneNode[2])
        whiteSpace = " " * (30 - len(name_childSpeciesNode_of_orthogroup))
        if targetGeneNode[2] == "NoGeneNode":
            resultLine = name_childSpeciesNode_of_orthogroup + "  " + whiteSpace + "NoGeneNode" + "  " + "NONE"
        else:
            resultLine = name_childSpeciesNode_of_orthogroup + "  " + whiteSpace + make_bsvalue_from_nodeLavel_NHXstyle(targetGeneNode[2]) + "  " + make_duplicationStatus_from_nodeLavel_NHXstyle(targetGeneNode[2])
        fSum.write(resultLine + "\n")
    fSum.write("\n")


    fSum.write(">SisterGeneGroups\n")
    #for targetSpeciesNode in speciesNodes_including_querySpecies:
    #for childSpeciesNode_of_orthogroup in childSpeciesNodes_orthogorup:
    #    if not name_querySpecies in childSpeciesNode_of_orthogroup[1]:
    #        continue
    for childSpeciesNode_of_orthogroup in childSpeciesNodes_orthogroup_including_querySpecies:
        name_childSpeciesNode_of_orthogroup = make_nodeName_from_nodeLavel_NHXstyle(childSpeciesNode_of_orthogroup[2])
        targetGeneNode = identify_targetGeneNode(allGeneNodesSR, name_childSpeciesNode_of_orthogroup, "SisterGeneGroups", topHitName_1stQuery)
        whiteSpace = " " * (30 - len(name_childSpeciesNode_of_orthogroup))
        if targetGeneNode[2] == "NoGeneNode":
            resultLine = name_childSpeciesNode_of_orthogroup + "  " + whiteSpace + "NONE   NoGeneNode"
        else:
            parentNode_queryGeneGroup, sisterGeneGroup = identify_sisterGeneNode(allGeneNodesSR, targetGeneNode)
            if sisterGeneGroup[2] == "NoSisterNode":
                resultLine = name_childSpeciesNode_of_orthogroup + "  " + whiteSpace + "NONE   NoSisterNode"
            else:
                nodeName_sisterGeneGroup = check_sisterGroupName_included_in_ancestralSpeciesNodeNames(allGeneNodesSR, childSpeciesNode_of_orthogroup, sisterGeneGroup)
                resultLine = name_childSpeciesNode_of_orthogroup + "  " + whiteSpace + make_bsvalue_from_nodeLavel_NHXstyle(parentNode_queryGeneGroup[2]) + "   " + nodeName_sisterGeneGroup
        fSum.write(resultLine + "\n")
    fSum.write("\n")
    ####

    #fSum.write(">BootstrapValue_sisterGeneGroup\n")
    #fSum.write(make_bsvalue_from_nodeLavel_NHXstyle(sisterGeneGroup[2]) + "\n\n")
    #fSum.write(">Members_sisterGeneGroup\n")
    #for member in sisterGeneGroup[1]:
    #    fSum.write(member + "\n")
    #fSum.write("\n")

    fSum.write(">Number_of_duplicatedNode\n")
    recs_duplications_for_speciesNodes = count_duplications_for_speciesNodes(allGeneNodesSR, topHitName_1stQuery)
    for nodeName, numDup in recs_duplications_for_speciesNodes.items():
        whiteSpace =  " " * (30 - len(nodeName)) 
        fSum.write(nodeName + "  " + whiteSpace + str(numDup) + "\n")
    fSum.write("\n")

    fSum.write(">AnalysisTime\n")
    elapsed_time = round((time.time() - startTime),1)
    fSum.write (str(elapsed_time) + " seconds")
    fSum.write("\n")

    fSum.close()


def makeSummary(outfile_summary):
    #print("#### makeSummary ####")
    recs_cds_assigned_by_ID = readFasta_dict(eachDirAddress, "000_cds_assigned_by_ID.txt")
    #cDNAfn = ""
    #cDNAfn = readFasta_dict(eachDirAddress, "000_cds_assigned_by_ID.txt")
    recAAfn = readFasta_dict(eachDirAddress, "000_translated_cds_assigned_by_ID.txt")

    rec044_unambSiteRate = ""
    fTMP = open(eachDirAddress + "040_mafOutAA.txt")
    fMafOut = list(fTMP)
    fTMP.close()
    if fMafOut:
        rec044_unambSiteRate = readFasta_dict(eachDirAddress, "044_aligned_site_rate.txt")
    
    fs = open(eachDirAddress + outfile_summary, "w")

    fs.write("################ Results: 1st analysis ################\n\n")

    fs.write(">QuerySequence\n")
    lines_query = make_lines_hit_query(recs_cds_assigned_by_ID)
    for line in lines_query:
        fs.write(line)
    fs.write("\n")

    fs.write(">Number_of_blastHits\n")
    rec_blastHits = count_blastHits(infile = "010_blastRes.txt")
    #print("rec_blastHits", rec_blastHits)
    #exit()
    rec_blastHits = whiteSpaceAdd(rec_blastHits)
    for name, num in rec_blastHits.items():
        fs.write(name + str(num) + "\n")
    fs.write("\n")

    speciesTree = ""
    f1stTree = open(eachDirAddress + "000_speciesTree.txt")
    speciesTree = "".join(list(f1stTree))
    speciesTree = re.sub("[ \n]", "", speciesTree)
    f1stTree.close()

    rearranged_1st_gene_tree_NHX = ""
    if os.path.exists(eachDirAddress + "085_NJBS1st.txt.rearrange.0"):
        f1stTree = open(eachDirAddress + "085_NJBS1st.txt.rearrange.0")
        rearranged_1st_gene_tree_NHX = list(f1stTree)[0]
    f1stTree.close()


    if not fMafOut:
        fs.write(">BS_of_orthogroupBasalNode\n")
        fs.write("No mafft out.\n")
        fs.write("\n")
    elif not rearranged_1st_gene_tree_NHX:
        fs.write(">BS_of_orthogroupBasalNode\n")
        fs.write("1st tree not estimated.\n")
        fs.write("\n")
    else:
        orthogroup = identify_orthogroup(rearranged_1st_gene_tree_NHX)
        if orthogroup[0] == 0:
            fs.write(">BS_of_orthogroupBasalNode\n")
            #if orthogroup[2] == "noOrthogroup_noKeynode":
            #    fs.write(orthogroup[2] + ": " + keyNode +"\n")
            #elif orthogroup[2] == "noOrthogroup_noQuerySequence":
            #    querySeq1 = lines_query[0]
            #    querySeq1 = re.sub(" <= .*$", "", querySeq1)
            #    #fs.write("Query, " + querySeq1 + ", was removed by the ShortSequence_threshold option\n")
            #    fs.write(orthogroup[2] + ": " + querySeq1 +"\n")
            #else:
            #    print("Error. Check orthogroup = identify_orthogroup(rearranged_1st_gene_tree_NHX).")
            #    exit()
            fs.write(orthogroup[2] + "\n")
            fs.write("\n")
        else:
            fs.write(">BS_of_orthogroupBasalNode\n")
            if len(orthogroup[1]) < 4:
                fs.write("Less than 4 orthogroup members.\n")
            else:
                fs.write(make_bsvalue_from_nodeLavel_NHXstyle(orthogroup[2]) + "\n\n")
            fs.write("\n")
    
            fs.write(">Orthogroup\n")
            sorted_members_focalGeneNode = sorted(orthogroup[1])
            for leaf in sorted_members_focalGeneNode:
                fs.write(leaf + "\n")
            fs.write("\n")
    
            fs.write(">GeneNumber_of_orthogroup\n")
            speciesNames_in_orthogroup = collect_speciesNames_in_orthogroup()
            #print("speciesNames_in_orthogroup", speciesNames_in_orthogroup)
            #exit()
            recs_geneNumber_of_orthogroup = OrderedDict()
            for species in speciesNames_in_orthogroup:
                count_hit = 0
                for leaf_gene in orthogroup[1]:
                    if re.search(species + "_", leaf_gene):
                        count_hit += 1
                recs_geneNumber_of_orthogroup[species] = count_hit
                #whiteSpece = " " * (30 - len(species))
                #fs.write(species + whiteSpece + str(count_hit) + "\n")
            #fs.write("\n")
            #print("recs_geneNumber_of_orthogroup", recs_geneNumber_of_orthogroup)
            recs_geneNumber_of_orthogroup = whiteSpaceAdd(recs_geneNumber_of_orthogroup)
            for species, count_hit in recs_geneNumber_of_orthogroup.items():
                fs.write(species + str(count_hit) + "\n")
            fs.write("\n")
    
            fs.write(">Rooting\n")
            rootGeneLeaves = selectRootSp4secondTreeSearch()
            for leaf in rootGeneLeaves:
                fs.write(leaf + "\n")
            fs.write("\n")
    
            #fs.write(">BootstrapValue_parentNode\n")
            #fs.write(make_bsvalue_from_nodeLavel_NHXstyle(bsValue_parentNode) + "\n\n")
    
        fs.write("\n")
    
        fs.write(">1st_rearranged_gene_tree_newick\n")
        rearrangedTreeNewick = change_nhx_to_newick_with_NHXnodeName(rearranged_1st_gene_tree_NHX)
        fs.write(rearrangedTreeNewick)
        fs.write("\n")
    
        fs.write(">1st_rearranged_gene_tree_NHX\n")
        fs.write(rearranged_1st_gene_tree_NHX)
        fs.write("\n")
    
        fs.write(">1st_gene_tree_newick\n")
        f1stTree = open(eachDirAddress + "085_NJBS1st.txt")
        fs.write(list(f1stTree)[0])
        f1stTree.close()
        fs.write("\n")
    
        if rec044_unambSiteRate:
            rec044_unambSiteRate = whiteSpaceAdd(rec044_unambSiteRate)
            fs.write(">Aligned-ShortSequence_threshold evaluation\n")
            for name, siteRate in rec044_unambSiteRate.items():
                nameRate = name + "  " + str(siteRate)
                if float(siteRate) > float(Aligned_site_rate):
                    fs.write(nameRate + "\n")
                else:
                    fs.write("== Removed ==> " + nameRate + "\n")
            fs.write("\n")


    fs.write("\n################ Settings ################\n\n")

    fs.write(">NumberAssigned_querySequence\n")
    for name, seq in recs_cds_assigned_by_ID.items():
        name = re.sub("[\n\r]", "", name)
        fs.write(name[1:] + "\n")
        fs.write(seq + "\n")
    fs.write("\n")

    fs.write(">SpeciesTree\n")
    fs.write(SpeciesTree + "\n")
    fs.write("\n")

    #fs.write(">SpeciesTree\n")
    #fs.write(speciesTree)
    #fs.write("\n\n")

    #fs.write(">KeyNode\n")
    #fs.write(keyNode)
    #fs.write("\n\n")

    if dataset == "Exclude3rd":
        fs.write(">SubstitutionModel\n" + "F84 (Tamura and Nei 1993) + gamma\n\n")
    elif dataset == "Include3rd":
        fs.write(">SubstitutionModel\n" + "TN93 (Tamura and Nei 1993) + gamma\n\n")
    else:
        fs.write(">SubstitutionModel\n" + "WAG (Whelan and Goldman 2001) + gamma\n\n")

    fs.write(">ShortSequence_threshold\n" + str(Aligned_site_rate) + "\n\n")

    fs.write(">Dataset\n" + dataset +  "\n")
    fs.write("\n")

    fs.write(">Rearrangement_BS_value_threshold\n")
    fs.write(str(BSthreshold) + "\n")
    fs.write("\n")


    fs.write(">TaxonSampling_color\n")
    for spNameTMP in taxonSamplingList:
        fs.write(spNameTMP + "\n")
    fs.write("\n")
    fs.close()
    


### makeSummary End
###############################################################


###############################################################
### Blast Start
def uniqueList(list_2d):
    recs_name_uniq  = OrderedDict()
    for ele in list_2d:
        nameLine = ele[0]
        identity = ele[1]
        if nameLine not in recs_name_uniq.keys():
            recs_name_uniq[nameLine] = identity
    
    return recs_name_uniq

def select_blastHitUnique(blastResFileFN):
    f = open(eachDirAddress + blastResFileFN)
    blastResAllLines = list(f)
    f.close()
    nameLines_all = []
    nameLine_tmp = ""
    flag = 0
    for line in blastResAllLines:
        line = line.rstrip("\n")
        if line.startswith(" Score") or line.startswith("Length="):
            continue
        if line.startswith(" Identities"):
            nameLine_tmp = re.sub("^> ", ">", nameLine_tmp)
            match = re.search("^ (Identities = [^,]+),", line)
            identity_line = match.group(1)
            nameLines_all.append([nameLine_tmp, identity_line])
            nameLine_tmp = ""
            flag = 0
        if flag == 1:
            nameLine_tmp += line
        if line.startswith(">"):
            flag = 1
            nameLine_tmp += line
    return uniqueList(nameLines_all)


def change_prohibitedExpression_in_nameLine(nameLine_FN):
    #nameLine_FN = re.sub("\+\+\-", "-", nameLine_FN)
    nameLine_FN = re.sub("\+\+", "+", nameLine_FN)
    nameLine_FN = re.sub("-+", "-", nameLine_FN)
    return nameLine_FN

def shorten_nameLine(nameLineTMP):
    #print("### shorten_nameLine ###")
    nameLine = re.sub("(>.{60}).*", r"\1", nameLineTMP)
    #print("nameLine", nameLine)
    return nameLine

def hitRecPicker():
    #print("#### hitRecPicker ####")
    blastResOut = open(eachDirAddress + "010_blastRes.txt", "w")
    AAout = open(eachDirAddress + "030_retrievedAAfas.txt",  "w")
    CDNAout = ""
    CDNAout = open(eachDirAddress + "030_retrievedDNAfas.txt", "w")

    #'''
    num_BlastpHits = 0
    for dbline in dbLines:

        recdbAA = readFasta_dict(dbAddress, dbline[1])
        recdbDNA = ""
        recdbDNA = readFasta_dict(dbAddress, dbline[2])

        blastResFile   = "005_vs" + dbline[0][:-1] + ".txt"
        #print("blastResFile:", blastResFile)
        recs_nameLine_identity_unique = select_blastHitUnique(blastResFile)
        
        for nameLine_blasthit in reversed(recs_nameLine_identity_unique.keys()):
            if not nameLine_blasthit.startswith(">"):
                continue
            
            num_BlastpHits += 1

            nameLine_blasthit = nameLine_blasthit.rstrip("\n")
            nameLine_blasthitTMP = re.sub(">", ">" + dbline[0], nameLine_blasthit)
            blastResOut.write(nameLine_blasthitTMP + "\n")

            counter_DBNLINEnum = 0
            for name_recdbAA, seq_recdbAA in recdbAA.items():
                if name_recdbAA == nameLine_blasthit:
                    break
                counter_DBNLINEnum += 1
            #match = re.search("DBNLINE\|(\d+)\|", nameLine_blasthit)
            #DBNLINEnum = match.group(1)
            DBNLINEnum = counter_DBNLINEnum
            nameLine = list(recdbAA.keys())[int(DBNLINEnum)]
            nameLine = change_prohibitedExpression_in_nameLine(nameLine)
            nameLine = re.sub(">", ">" + dbline[0][:-1] + "_", nameLine)
            #print("nameLine",nameLine)
            
            #nameLine = re.sub("(>.{60}).*", r"\1", nameLine)
            nameLine = shorten_nameLine(nameLine)
            
            #print("nameLine",nameLine)
            #print("")
            AAout.write(nameLine + "\n")
            AAout.write(list(recdbAA.values())[int(DBNLINEnum)] + "\n")
            #CDNAout.write(list(recdbDNA.keys())[int(DBNLINEnum)]   + "\n") # modified 20241118
            CDNAout.write(nameLine + "\n") # modified 20241118
            CDNAout.write(list(recdbDNA.values())[int(DBNLINEnum)] + "\n")
            
            #print()
    blastResOut.close()

    if num_BlastpHits < 4:
        result = "Less than 4 blast hits."
        error_makeSummary(result)
        #print (result + " <br>")
        if Switch_deleteIntermediateFiles == "L":
            error_resHtmlMaker(result)
        if Switch_deleteIntermediateFiles == "D":
            deleteFiles()
        exit()

    AAout.close()
    CDNAout.close()


def make_spNames_IndbLines():
    spNames_IndbLines = []
    for dbLine in dbLines:
        spNamePrefix_dbLine = dbLine[0]
        spNamePrefix_dbLine = re.sub("_$", "", spNamePrefix_dbLine)
        spNames_IndbLines.append(spNamePrefix_dbLine)
    return spNames_IndbLines


def blastpSearch():
    DirAafileName = eachDirAddress + "000_translated_cds_assigned_by_ID.txt"

    for dbline in dbLines:
        dbAAfile = dbAddress + dbline[1]
        outFile = eachDirAddress + "005_vs" + dbline[0][:-1] + ".txt"
        comLine = "tools/blastp -query {0} -evalue {1} -num_alignments {2}  -num_descriptions {3}  -db {4} -out {5}"\
                  .format(DirAafileName, blastEvalue, Number_of_hits_to_report_per_genome, Number_of_hits_to_report_per_genome, dbAAfile, outFile)
        #print("comLine", comLine)
        subprocess.call(comLine, shell=True)


def make_lines_hit_query(recs_cds_assigned_by_ID):
    topHits = topHitPicker()

    topHitValue0s = [x[0] for x in topHits.values()]
    length_longestName = len(max(topHitValue0s, key = len))

    lines_hit_query_FN = []
    for i in range(len(topHits)):
        uploadedSeqName = list(recs_cds_assigned_by_ID.keys())[i][1:]
        uploadedSeqName = re.sub("[\n\r]", "", uploadedSeqName)
        InfoIdentity = list(topHits.values())[i][1]
        if not InfoIdentity:
            InfoIdentity = "Name"
        lines_hit_query_FN.append(list(topHits.values())[i][0][1:] \
                 + " " * (length_longestName - len(list(topHits.values())[i][0][1:])) \
                 + " <= " \
                 + "[" + InfoIdentity + "] "\
                 #+ " " * (30 - len(InfoIdentity)) \
                 + ": "+ uploadedSeqName \
                 + "\n")
    return lines_hit_query_FN


def topHitPicker():
    topHits = OrderedDict()
    recs_cds_assigned_by_ID = readFasta_dict(eachDirAddress, "000_cds_assigned_by_ID.txt")
    for nameLine_cds_assigned_by_ID in recs_cds_assigned_by_ID.keys():
        match = re.search(">([^_]+)_", nameLine_cds_assigned_by_ID)
        querySpecies_SR = match.group(1)
        #recs_blastnRes = read_blastnRes2(eachDirAddress, "005_vs" +querySpecies_SR + ".txt")
        recs_blastnRes = select_blastHitUnique("005_vs" + querySpecies_SR + ".txt")
        blastTopHitNameLine = ""
        blastTopHitIdentity = ""
        for nameLine_blasthit, identity_blasthit in recs_blastnRes.items():
            blastTopHitNameLine = re.sub(" .*$", "", nameLine_blasthit)
            #print("blastTopHitNameLine", blastTopHitNameLine)
            blastTopHitNameLine = change_prohibitedExpression_in_nameLine(blastTopHitNameLine)
            #print("blastTopHitNameLine", blastTopHitNameLine)
            print("")
            blastTopHitNameLine = re.sub(">", ">" + querySpecies_SR + "_", blastTopHitNameLine)
            blastTopHitIdentity = identity_blasthit
            break
        blastTopHitNameLine = shorten_nameLine(blastTopHitNameLine)
        topHits[nameLine_cds_assigned_by_ID] = [blastTopHitNameLine, blastTopHitIdentity.rstrip("\n")]
    return topHits


### Blast End
###############################################################


###############################################################
### delete_sequences_with_alignedSiteRate Start
def compare_query2other_aliSiteRate(querySeq, otherSeq):
    # This function is not used and replaced with compare_query2other_aliSiteRate_nogap.
    count = 0
    for i in range(len(querySeq)):
        if re.match("\w", querySeq[i]) and re.match("\w", otherSeq[i]):
            count += 1
        elif re.match("-", querySeq[i]) and re.match("\w", otherSeq[i]):
            count += 1
        else:
            continue
    return round(float(count)/len(querySeq),3)

def compare_query2other_aliSiteRate_nogap(querySeq, otherSeq):
    count = 0
    for i in range(len(querySeq)):
        if re.match("-", querySeq[i]):
            continue
        if re.match("\w", otherSeq[i]):
            count += 1
        else:
            continue
    querySeq_noGap = re.sub("-", "", querySeq)
    #print("querySeq", querySeq)
    #print("querySeq_noGap", querySeq_noGap)
    #exit()
    return round(float(count)/len(querySeq_noGap),3)


def calculate_nonGapSiteRate(trimledAAfile):
    #print("#### calculate_nonGapSiteRate ####")
    recs_nonGapSiteRateFN = OrderedDict()
    recAA = readFasta_dict(eachDirAddress, trimledAAfile)
    #print("trimledAAfile", trimledAAfile)
    #exit()
    #for name, seq in recAA.items():
    #    print("name", name)
    #exit(9)

    querySeq = list(recAA.values())[-1]
    queryName = list(recAA.keys())[-1]
    #print("queryName", queryName)
    #exit()
    for otherName, otherSeq in recAA.items():
        protID = re.sub(" .*$","", otherName)
        #aaRate = compare_query2other_aliSiteRate(querySeq, otherSeq)
        aaRate = compare_query2other_aliSiteRate_nogap(querySeq, otherSeq)
        #print("protID", protID)
        #print("aaRate", aaRate)
        #print("querySeq", querySeq)
        #print("otherSeq", otherSeq)
        recs_nonGapSiteRateFN[protID] = aaRate
    #exit()
    return recs_nonGapSiteRateFN


def delete_sequences_with_alignedSiteRate():
    #print("### Now fixing: delete_sequences_with_alignedSiteRate 1 ###")

    recs_nonGapSiteRate = calculate_nonGapSiteRate("042_AA.fas.trm")
    #exit()
    
    #for name, rate in recs_nonGapSiteRate.items():
    #    print("name1", name)
    #    print("rate1", rate)
    #print("### delete_sequences_with_alignedSiteRate 2 ###")

    file_overRateAA = open(eachDirAddress + "044_overRateAA.fas", "w")
    file_overRateDNA = ""
    file_overRateDNA = open(eachDirAddress + "044_overRateDNA.fas", "w")
    file_unambSiteRatefile = open(eachDirAddress + "044_aligned_site_rate.txt", "w")

    recAA = readFasta_dict(eachDirAddress, "030_retrievedAAfas.txt")
    recDNA = ""
    recDNA = readFasta_dict(eachDirAddress, "030_retrievedDNAfas.txt")

    #print("len(recs_nonGapSiteRate)", len(recs_nonGapSiteRate))
    #print("len(recAA)",len(recAA))
    #exit()
    if len(recs_nonGapSiteRate) == len(recAA):
        for i in range (len(recs_nonGapSiteRate)):
            name = list(recs_nonGapSiteRate.keys())[i]
            rate = list(recs_nonGapSiteRate.values())[i]
            #print("name1", name)
            #print("rate1", rate)
            #print("Aligned_site_rate", float(Aligned_site_rate))
            if rate > float(Aligned_site_rate):
                print(rate, Aligned_site_rate)
                #print("(list(recAA.keys())[i]", i, list(recAA.keys())[i]) # modified 20241118
                file_overRateAA.write  (list(recAA.keys())[i]    + "\n")
                file_overRateAA.write  (list(recAA.values())[i]  + "\n")
     
                #print("(list(recDNA.keys())[i]", i, list(recDNA.keys())[i])
                file_overRateDNA.write(list(recDNA.keys())[i]   + "\n") # sequence names should be unique 
                file_overRateDNA.write(list(recDNA.values())[i] + "\n")
    
            file_unambSiteRatefile.write(name + "\n" + str(rate) + "\n")
            #print("")
    else:
        ### 042_AA.fas.trm  : Only gap sequence is moved by TRIMAL1.4
        for name_siterate, rate_siterate in recs_nonGapSiteRate.items():
            #print("name_siterate", name_siterate)
            match_nsr = re.search("^[^_]+_([^_]+)_", name_siterate)
            name_siterate_geneID = match_nsr.group(1)
            #print("name_siterate_geneID", name_siterate_geneID)
            #print("rate_siterate", rate_siterate)
            if float(rate_siterate) > float(Aligned_site_rate):
                
                for name_AA, seq_AA in recAA.items():
                    #print("name_AA", name_AA)
                    if re.search("_" + name_siterate_geneID + "_", name_AA):
                        file_overRateAA.write(name_AA + "\n")
                        file_overRateAA.write(seq_AA + "\n")
                        break

                for name_DNA, seq_DNA in recDNA.items():
                    #print("name_DNA", name_DNA)
                    if re.search(">" + name_siterate_geneID + " ", name_DNA):
                        file_overRateDNA.write(name_DNA + "\n")
                        file_overRateDNA.write(seq_DNA + "\n")
                        break
    
            file_unambSiteRatefile.write(name_siterate + "\n" + str(rate_siterate) + "\n")
            #print("")
    file_overRateAA.close()
    file_overRateDNA.close()
    file_unambSiteRatefile.close()
    #exit()

### delete_sequences_with_alignedSiteRate End
###############################################################


###############################################################
### alignmentFile_PhyAnal Start
def orderedDict2FasFile(recs, outfile):
    out           = open(eachDirAddress + outfile, "w")
    for name,value in recs.items():
        out.write(name + "\n")
        out.write(value + "\n")
    out.close()


def orderedDict2phyFile(recs, outfile):
    secLength     = len(sorted(recs.values())[0])
    spSeqSizeLine = str(len(recs)) + " " + str(secLength)

    recs = whiteSpaceAdd(recs)
    out           = open(eachDirAddress + outfile, "w")
    out.write(spSeqSizeLine + "\n")
    for name,value in recs.items():
        out.write(name + value + "\n")
    out.close()


def read_TrimalHTMLout_dict_v12(trimAAresult):
    f = open(eachDirAddress + trimAAresult)
    lines = list(f)
    f.close()

    recs_trimalHTMLout  = OrderedDict()
    for line in lines:
        if re.search("    <span class=sel>Selected Sequences",line):
            continue
        if line.startswith("    <span class=sel>"):
            line = line.rstrip("\n")
            match = re.search("<span class=sel>([^<]+)</span> +([^ ].*)$", line)
            name     = match.group(1)
            sequence = match.group(2)
            if not name in recs_trimalHTMLout.keys():
                recs_trimalHTMLout[name] = ""
            recs_trimalHTMLout[name] += sequence

    return recs_trimalHTMLout

def read_TrimalHTMLout_dict_v141(trimAAresult):
    f = open(eachDirAddress + trimAAresult)
    lines = list(f)
    f.close()

    trimalMarkedSites = ""
    for line in lines:
        if line.startswith("    Selected Cols:"):
            line = line.rstrip("\n")
            sequence = re.sub(" +Selected Cols: +", "",line)
            trimalMarkedSites += sequence

    return trimalMarkedSites


def trimaledv12_FileMakerDNA(fastaFile, trimAAresult, outfile):
    recs = readFasta_dict(eachDirAddress, fastaFile)
    
    recs_trimalHTMLout = read_TrimalHTMLout_dict_v12(trimAAresult)
    #print("recs_trimalHTMLout", recs_trimalHTMLout)
    #exit()
    
    trimalMarkedSites = list(recs_trimalHTMLout.values())[-1]
    #print("trimalMarkedSites", trimalMarkedSites)
    #exit()
    trimalMarkedSites = re.sub("<span class=sel>.</span>", "#", trimalMarkedSites)
    #print("trimalMarkedSites", trimalMarkedSites)
    #exit()

    recsTrimed = OrderedDict()
    for name,value in recs.items():
        sequence = ""
        for i in range(len(trimalMarkedSites)):
            if trimalMarkedSites[i] == "#":
                #out.write(trimalMarkedSites[i])
                sequence += value[i*3] + value[i*3+1] + value[i*3+2]
        recsTrimed[name] = sequence
    orderedDict2phyFile(recsTrimed, outfile = outfile)


def trimaledv141_FileMakerDNA(fastaFile, trimAAresult, outfile):
    recs = readFasta_dict(eachDirAddress, fastaFile)
    #print("fastaFile", fastaFile)
    #exit()
    
    trimalMarkedSitesTMP = read_TrimalHTMLout_dict_v141(trimAAresult)
    #print("trimalMarkedSitesTMP", trimalMarkedSitesTMP)
    #exit()
    trimalMarkedSites = re.sub("<span class=nsel> </span>", "-", trimalMarkedSitesTMP)
    trimalMarkedSites = re.sub("<span class=sel> </span>", "#", trimalMarkedSites)
    #print("trimalMarkedSites", trimalMarkedSites)
    #exit()

    recsTrimed = OrderedDict()
    for name,value in recs.items():
        sequence = ""
        for i in range(len(trimalMarkedSites)):
            if trimalMarkedSites[i] == "#":
                #out.write(trimalMarkedSites[i])
                sequence += value[i*3] + value[i*3+1] + value[i*3+2]
        recsTrimed[name] = sequence
    orderedDict2phyFile(recsTrimed, outfile = outfile)


### alignmentFile_PhyAnal End
###############################################################


###############################################################
### alignmentFile_HTML
def outGroupSelect(phyFileName):
    recSeqFN = readPhy_dict(phyFileName)
    outgroupTMP = list(recSeqFN.keys())[0]
    return outgroupTMP[1:]


def reorderSeqByTree(recsFN, treeFileName):
    leaves = []
    leaves = collect_leaves_InOrderFrom_bothNHXnewick(treeFileName)
    seqDictFN  = OrderedDict()
    for leaf in reversed(leaves):
        Lleaf = ">" + leaf
        seqDictFN[Lleaf] = recsFN[Lleaf]
    return seqDictFN


def delete_nameSpaceSeqBp(recsFN):
    recsFN2        = OrderedDict()
    for name, seq in recsFN.items():
        if re.search(" \d+ bp$", name):
            name = re.sub(" \d+ bp$", "", name)
        recsFN2[name] = seq
    return recsFN2


def nameChange_whiteLaterDelete(recsFN):
    recsFN2        = OrderedDict()
    for name, seq in recsFN.items():
        name = re.sub(" .*$", "", name)
        recsFN2[name] = seq
    return recsFN2


def gapDelete(recsFN):
    recsFN2        = OrderedDict()
    for name, seq in recsFN.items():
        seq = re.sub("-", "", seq)
        recsFN2[name] = seq
    return recsFN2


def whiteSpaceAdd(recsFN1):
    longestName = max(recsFN1.keys(), key = len)
    longestName = re.sub("<[^>]+>", "", longestName)
    longestNameLen = len(longestName)
    recsFN2        = OrderedDict()
    for name,value in recsFN1.items():
        name = re.sub("^>", "", name)
        nameTMP = re.sub("<[^>]+>", "", name)
        nameWhiteSpace = name + " " * (longestNameLen - len(nameTMP) + 2)
        recsFN2[nameWhiteSpace] = value
    return recsFN2


def reorderDeleteGap_Fas2FasByTree(fastaFileName, tree4leafOrder, outPhyFileName):
    recs = readFasta_dict(eachDirAddress,fastaFileName) 
    recs = nameChange_whiteLaterDelete(recs)
    recs = gapDelete(recs)
    recs = reorderSeqByTree(recs, tree4leafOrder)
    orderedDict2FasFile(recs, outfile = outPhyFileName)


def phy2fastmePhy(phyFileName, outFastmePhyFileName):
    recsFN = readPhy_dict(phyFileName)
    secLength     = len(sorted(recsFN.values())[0])
    spSeqSizeLine = str(len(recsFN)) + " " + str(secLength)

    recsFN = whiteSpaceAdd(recsFN)
    out  = open(eachDirAddress + outFastmePhyFileName, "w")
    out.write(spSeqSizeLine + "\n")
    for name,value in reversed(recsFN.items()):
        out.write(name + value + "\n")
    out.close()


def fas2phy(fastaFileName, outPhyFileName):
    recs = readFasta_dict(eachDirAddress,fastaFileName)
    recs = delete_nameSpaceSeqBp(recs)
    orderedDict2phyFile(recs, outfile = outPhyFileName)


def reorderFas2PhyByTree(fastaFileName, tree4leafOrder, outPhyFileName):
    recs = readFasta_dict(eachDirAddress,fastaFileName) 
    recs = delete_nameSpaceSeqBp(recs)
    for name, seq in recs.items():
        print(name)
        print("<br><br>")
    exit()
    recs = reorderSeqByTree(recs, tree4leafOrder)
    orderedDict2phyFile(recs, outfile = outPhyFileName)


def codonSepalate(recs):
    recs1 = OrderedDict()
    recs2 = OrderedDict()
    recs3 = OrderedDict()
    for name, sec in recs.items():
        recs1[name] = ""
        recs2[name] = ""
        recs3[name] = ""
        for i in range(len(sec)):
            if   i%3 == 0:
               recs1[name] += sec[i]
            elif i%3 == 1:
               recs2[name] += sec[i]
            else:
               recs3[name] += sec[i]
    return recs1, recs2, recs3


def phyCodonToBlock(phyFileName, blockNum, outfile):
    recs = readPhy_dict(phyFileName)
    recs1, recs2, recs3 = codonSepalate(recs)
    recsS1 = OrderedDict()
    if blockNum == 3:
        for name in recs1.keys():
            recsS1[name] = recs1[name] + recs2[name] + recs3[name]
            
    if blockNum == 2:
        for name in recs1.keys():
            recsS1[name] = recs1[name] + recs2[name]

    recsS2 = whiteSpaceAdd(recsS1)

    out = open(eachDirAddress + outfile, "w")
    out.write(str(len(recsS2)) + " " + str(len(list(recsS2.values())[0])) + "\n")
    for name,sec in recsS2.items():
        out.write(name + sec + "\n")
    out.close()


def identify_rootGeneLeaf_4_2ndAnalysis_when_noParantNode(allGeneNode_SR, focalGeneNode_SR):
    daughterNode_1st, daughterNode_2nd = identify_daughterNodes(allGeneNode_SR, focalGeneNode_SR)
    if daughterNode_1st[1] <= daughterNode_2nd[1]:
        rootGeneLeaf = list(daughterNode_1st[1])[0]
    else:
        rootGeneLeaf = list(daughterNode_2nd[1])[0]
    return [rootGeneLeaf]


def selectRootSp4secondTreeSearch():
    f1stTree = open(eachDirAddress + "085_NJBS1st.txt.rearrange.0")
    treeNHX = list(f1stTree)[0]
    f1stTree.close()
    allGeneNodes = collect_nodes_from_NHX(treeNHX)
    orthogroup = identify_orthogroup(treeNHX)
    list_rootSpecies = []
    if not (orthogroup):
        list_rootSpecies = identify_rootGeneLeaf_4_2ndAnalysis_when_noParantNode(allGeneNodes, orthogroup)
    else:
        #sisterGeneGroups_gettingDeeper = collect_sisterGroups(allGeneNodes, focalGeneGroup)   ##### Select rooting within orthogroup
        sisterGeneGroups_gettingDeeper = collect_sisterGroups(allGeneNodes, orthogroup)   ##### Select rooting outside of orthogroup

        candidates_rootSpecies = []
        for sisterGeneGroup in sisterGeneGroups_gettingDeeper:
            candidates_rootSpecies += list(sisterGeneGroup[1])

        if not candidates_rootSpecies:
            daughterNode_1st, daughterNode_2nd = identify_daughterNodes(allGeneNodes, orthogroup)
            topHits = topHitPicker()
            topHitName_1stQuery = list(topHits.values())[0][0]
            focalGeneGroup = ""
            rootingGeneGroup = ""
            if topHitName_1stQuery[1:] in daughterNode_1st[1]:
                rootingGeneGroup = daughterNode_2nd
                focalGeneGroup = daughterNode_1st
            else:
                rootingGeneGroup = daughterNode_1st
                focalGeneGroup = daughterNode_2nd
            candidates_rootSpecies = list(rootingGeneGroup[1])

        for i in range(0, int(num_rootSequences)):
            if i < len(candidates_rootSpecies):
                list_rootSpecies.append(candidates_rootSpecies[i])

    list_rootSpecies.reverse()
    return list_rootSpecies


def make_2ndanalysis_seqFile(rootLeaves_SR, outfile):
    recs_nucl = readFasta_dict(eachDirAddress, "054_p2nOutcDNAfas.txt")
    leaves_add_2_focalClade = []
    for rootLeaf in rootLeaves_SR:
       if not rootLeaf in resDict_1st[">Orthogroup"]:
            leaves_add_2_focalClade.append(rootLeaf)
    for leaf in resDict_1st[">Orthogroup"]:
        leaves_add_2_focalClade.append(leaf)

    out = open(eachDirAddress + "/" + outfile, "w")
    for name in leaves_add_2_focalClade:
        out.write(">" + name + "\n")
        seq = re.sub("-", "", recs_nucl[">" + name])
        out.write(seq + "\n")
    out.close()


def cDNAfas2noGapAAFasFile(cDNAfasFileName, outfile):
    recsFN = readFasta_dict(eachDirAddress, cDNAfasFileName)
    fa = open(eachDirAddress + "/" + outfile, "w")
    for name, seq in recsFN.items():
        fa.write(name + "\n")
        fa.write(translation(seq) + "\n")
    fa.close()


def make_raxmlPartitionFile(outPartFile):
    phyFile = open(eachDirAddress + "210_trimedBlockExc3rdPhy.txt", "r")
    phyLines = list(phyFile)
    phyFile.close() 
    spNumTMP, seqLength = re.split(" ",  phyLines[0])
    outgroup, seqTMP    = re.split(" +", phyLines[1])
    
    seqLength = int(seqLength)
    partFile = open(eachDirAddress + outPartFile, "w")
    partFile.write("DNA,gene1=1-"                                    + str(int(seqLength/2)) + "\n")
    partFile.write("DNA,gene2=" + str(int(seqLength/2)     + 1)     + "-" + str(int(seqLength))   + "\n")
    partFile.close()


def moveRAxMLfiles(outfile):
    line1 = "cp " + eachDirAddress + "RAxML_bipartitions.txt " + eachDirAddress + outfile
    #print(line1)
    subprocess.call(line1, shell=True)
    line2 = "rm " + eachDirAddress + "RAxML*"
    subprocess.call(line2, shell=True)

### alignmentFile_HTML End
###############################################################


###############################################################
### result_HTML Start
def error_resHtmlMaker(resultFN):
    topHits = topHitPicker()
    topHitName_1stQuery = list(topHits.values())[0]
    firstQueryTMP = topHitName_1stQuery[0][1:]
    firstQueryTMP = shorten_nameLine(firstQueryTMP)
    resHTMLlines1 = re.sub('FIRSTQUERY', firstQueryTMP, resHTMLlines_incomplete)
    if resultFN == "noOrthogroup_noQuerySequence":
        resultFN = resultFN + ": " + firstQueryTMP
    elif resultFN == "noOrthogroup_noKeynode":
        resultFN = resultFN + ": " + keyNode
    #print("resultFN111", resultFN)
    resHTMLlines1 = re.sub('BSVALUE_orthogroup_1STTREE', resultFN, resHTMLlines1)
    resHTMLlines1 = re.sub('EACHDIRADDRESS_', eachDirAddress, resHTMLlines1)

    #out = open(eachDirAddress + "300_resultsREA.html", "w")
    out = open(queryID + ".html", "w")
    out.write(resHTMLlines1)
    out.close()




def make_resHtml_link_form_outside():
    top = '''
    <!DOCTYPE html>
    <html>
    <head>
        <meta http-equiv="Content-Type" content="text/html">
            <title>TITLE</title>
                <style type="text/css">
                    .blackBG { background-color: #000000; color: white}
                    .redBG { background-color: #FF0000; color: white}
                </style>
        </head>
    <body>
    <pre>
    <span style="font-size: 130%;">'''
    
    bottom = '''
    </span></pre>
    </body>
    </html>'''

    out= open("draw_tree_" + queryID + ".html", "w")
    out.write(top + "\n")
    out.write('<a href="' + eachDirAddress + '/300_resultsREA.html" target="_blank">' + queryID + '</a>' + '\n')
    out.write(bottom + "\n")
    out.close()


def make_resHtml2(resHTMLlinesFN):
    
    #resDict_1stSummary = readRes_dict(eachDirAddress + "100_analysisSummary.txt")
    resDict_1st = readRes_dict(eachDirAddress + "/100_analysisSummary.txt")
    topHitName_1stQuery = resDict_1st[">QuerySequence"][0]
    topHitName_1stQuery = re.sub(" +.*", "", topHitName_1stQuery)
    firstQueryTMP = shorten_nameLine(topHitName_1stQuery)
    resHTMLlines1 = re.sub('FIRSTQUERY', topHitName_1stQuery, resHTMLlinesFN)

    #f1stTree = open(eachDirAddress + "085_NJBS1st.txt.rearrange.0")
    #treeNHX = list(f1stTree)[0]
    #f1stTree.close()

    #orthogroup = identify_orthogroup(treeNHX)
    
    ######### 1st tree
    #SisterNode_of_nodeIncludingQuery = ""
    #if ">SisterGeneGroup" in resDict_1st.keys():
    #    SisterNode_of_nodeIncludingQuery = resDict_1st[">SisterGeneGroup"][0]
    #else:
    #    SisterNode_of_nodeIncludingQuery = "Not_identified_sisterNode"
    #resHTMLlines1 = re.sub('SISTERNODE_VERATEBRATES', SisterNode_of_nodeIncludingQuery, resHTMLlines1)

    #BootstrapValue_SisterNode_of_nodeIncludingQuery = ""
    #if ">BootstrapValue_sisterGeneGroup" in resDict_1st.keys():
    #    BootstrapValue_SisterNode_of_nodeIncludingQuery = resDict_1st[">BootstrapValue_sisterGeneGroup"][0]
    #else:
    #    BootstrapValue_SisterNode_of_nodeIncludingQuery = "Not_identified_sisterNode"
    #resHTMLlines1 = re.sub('BS_2NDTREE_SISTERNODE', BootstrapValue_SisterNode_of_nodeIncludingQuery, resHTMLlines1)

    BootstrapValue_orthogroupBasalNode_1sttree = ""
    if ">BS_of_orthogroupBasalNode" in resDict_1st.keys():
        BootstrapValue_orthogroupBasalNode_1sttree = resDict_1st[">BS_of_orthogroupBasalNode"][0]
    else:
        BootstrapValue_orthogroupBasalNode_1sttree = "Not_identified_focalClade"
        #resHTMLlines1 = re.sub('<td colspan="2">Sister clade of vertebrate gene clade</td>', '<td>&nbsp;</td>', resHTMLlines1)
        resHTMLlines1 = re.sub('<td>Alignment: <a href="170_aln_prot.html" target="_blank">Amino acid</a>, <a href="190_aln_nucl.txt" target="_blank">Nucleotide</a></td>', '<td>Not analyzed.</td>', resHTMLlines1)
    resHTMLlines1 = re.sub('BSVALUE_orthogroup_1STTREE', BootstrapValue_orthogroupBasalNode_1sttree, resHTMLlines1)

    resHTMLlines1 = re.sub('EACHDIRADDRESS_', eachDirAddress, resHTMLlines1)


    #BootstrapValue_parentNode = ""
    #if ">BootstrapValue_parentNode" in resDict_1st.keys():
    #    BootstrapValue_parentNode = resDict_1st[">BootstrapValue_parentNode"][0]
    #else:
    #    BootstrapValue_parentNode = "No_orthogroup"
    #resHTMLlines1 = re.sub('BS_PARENT_1STTREE', BootstrapValue_parentNode, resHTMLlines1)

    ######## 2nd tree
    #BootstrapValue_nodeIncludingQuery = ""
    #if ">BootstrapValue_queryGeneGroup" in resDict_1st.keys():
    #    BootstrapValue_nodeIncludingQuery = resDict_1st[">BootstrapValue_queryGeneGroup"][0]
    #else:
    #    BootstrapValue_nodeIncludingQuery = "No_vertebrateNode"
    #resHTMLlines1 = re.sub('BS_2NDTREE_VERTEBRATENODE', BootstrapValue_nodeIncludingQuery, resHTMLlines1)

    #BootstrapValue_sisterNode_and_nodeIncludingQuery = ""
    #if ">BootstrapValue_sisterGeneGroup_vs_queryGeneGroup" in resDict_1st.keys():
    #    BootstrapValue_sisterNode_and_nodeIncludingQuery = resDict_1st[">BootstrapValue_sisterGeneGroup_vs_queryGeneGroup"][0]
    #else:
    #    BootstrapValue_sisterNode_and_nodeIncludingQuery = "No_orthogroup"
    #resHTMLlines1 = re.sub('BS_2NDTREE_PARENTNODE', BootstrapValue_sisterNode_and_nodeIncludingQuery, resHTMLlines1)

    #out = open(eachDirAddress + "300_resultsREA.html", "w")
    out = open(queryID + ".html", "w")
    out.write(resHTMLlines1)
    out.close()

### result_HTML End
###############################################################


###############################################################
### Others START

def deleteFiles():

    keywords_dirFiles = [
             "0",
             "1[1-9]",
             "2",
             "3",
            ]
    fileNames = os.listdir(path=eachDirAddress)
    fileNames_rm = []
    for fileName in fileNames:
        for keyword in keywords_dirFiles:
            if re.search("^" + keyword, fileName):
                fileNames_rm.append(fileName)
    for file in fileNames_rm:
        address_file = eachDirAddress + file
        line_rm = "rm " + address_file
        #print("line_rm", line_rm)
        subprocess.call(line_rm, shell=True)

### Others END


########### Data summarize after all gene tree estimated
def make_lines_atmarkSeparated(geneIDs_fn):

    speciesNames_in_orthogroup = collect_speciesNames_in_orthogroup()


    #querySpeciesNode = identify_speciesNode(name_querySpecies)
    #speciesNodes_including_querySpecies = collect_ancestralNodes(allNodes_speciesTree, querySpeciesNode)

    lines_FN = []
    for x in range(len(geneIDs_fn)):
        print(x+1, geneIDs_fn[x])
        
        line_FN = "Num" + "@" + str(x+1) + ","

        line_FN += "QueryGeneID" + "@" +geneIDs_fn[x] + ","

        address_eachDirectory = outdir + geneIDs_fn[x]
        if not os.path.exists(address_eachDirectory):
            print("Cannot find directories. Stopped:")
            print(address_eachDirectory)
            exit()

        address_100_resultsREA_html = outdir + geneIDs_fn[x] + "/100_analysisSummary.txt"
        if not os.path.exists(address_100_resultsREA_html):
            print("Cannot find 100_analysisSummary.txt. Stopped:")
            print(address_100_resultsREA_html)
            exit()

        address_100_analysisSummary_txt = outdir + "/" + geneIDs_fn[x] + "/100_analysisSummary.txt"
        if not os.path.exists(address_100_analysisSummary_txt):
            print("Cannot find 100_analysisSummary.txt. Stopped:")
            print(address_100_analysisSummary_txt)
            exit()

        seqDict = readRes_dict(address_100_analysisSummary_txt)

        #for name, val in seqDict.items():
        #    print(name)
        #    #print(val)
        #exit()
    
        line_FN += "QueryLength@"
        if ">NumberAssigned_querySequence" in  seqDict.keys():
            line_FN += str(len(seqDict[">NumberAssigned_querySequence"][1])) + ","
        else:
            line_FN += "NONE" + ","

        line_FN += "SpeciesWithGeneFunction@"
        if ">Orthogroup" in  seqDict.keys():
            flagTMP = 0
            for geneLeaf in seqDict[">Orthogroup"]:
                if re.search("^" + speciesWithGeneFunction + "_", geneLeaf):
                    flagTMP += 1
                    line_FN += geneLeaf + " "
            if flagTMP < 1:
                line_FN += "NONE,"
            else:
                line_FN += ","
            #exit()
        else:
            line_FN += "NONE" + ","
    
        line_FN += "BS_of_orthogroupBasalNode@"
        if ">BS_of_orthogroupBasalNode" in  seqDict.keys():
            line_FN += seqDict[">BS_of_orthogroupBasalNode"][0] + ","
        else:
            line_FN += "NONE" + ","

        line_FN += "2ndGeneTree@"
        if ">2nd_rearranged_gene_tree_newick" in  seqDict.keys():
            line_FN += "DONE" + ","
        else:
            line_FN += "NONE" + ","

        #line_FN += "BootstrapValue_sisterGeneGroup_vs_queryGeneGroup@"
        #if ">BootstrapValue_sisterGeneGroup_vs_queryGeneGroup" in  seqDict.keys():
        #    line_FN += seqDict[">BootstrapValue_sisterGeneGroup_vs_queryGeneGroup"][0] + ","
        #else:
        #    line_FN += "NONE" + ","
    
        #line_FN += "SisterGeneGroup@"
        #if ">SisterGeneGroup" in  seqDict.keys():
        #    print("seqDict[>SisterGeneGroup]", seqDict[">SisterGeneGroup"])
        #    exit()
        #    line_FN += seqDict[">SisterGeneGroup"][0] + ","
        #else:
        #    line_FN += "NONE" + ","

        
        if ">Number_of_blastHits" in  seqDict.keys():
            for node_num in seqDict[">Number_of_blastHits"]:
                match = re.search("^([^ ]+) +(\d+)$", node_num)
                node = match.group(1)
                num = match.group(2)
                line_FN += "BHnum_" + node + "@" + num + ","
        else:
            #for speciesName_in_orthogroup in speciesNames_in_orthogroup:
            for dbLine in dbLines:
                speciesName = dbLine[0]
                #print("speciesName_in_orthogroup", speciesName_in_orthogroup)
                #print("speciesName", speciesName[:-1])
                line_FN += "BHnum_" + speciesName[:-1] + "@NONE" + ","

        if ">GeneNumber_of_orthogroup" in  seqDict.keys():
            for node_num in seqDict[">GeneNumber_of_orthogroup"]:
                match = re.search("^([^ ]+) +(\d+)$", node_num)
                node = match.group(1)
                num = match.group(2)
                line_FN += "OGnum_" + node + "@" + num + ","
        else:
            for speciesName_in_orthogroup in speciesNames_in_orthogroup:
                #line_FN += "OGnum_" + speciesName_in_orthogroup + "@NONE" + ","
                #line_FN += "OGnum_" + speciesName_in_orthogroup + "@" + " " + ","
                line_FN += "OGnum_" + speciesName_in_orthogroup + "@NONE" + ","

        if ">MonophyleticGeneGroups" in  seqDict.keys():
            #print("seqDict[>Number_of_duplicatedNode]", seqDict[">Number_of_duplicatedNode"])
            for node_BS_sister in seqDict[">MonophyleticGeneGroups"]:
                #print("node_BS_sister", node_BS_sister, "|")
                match = re.search("^([^ ]+) +([^ ]+) +([^ ]+)$", node_BS_sister)   ##################
                #match = re.search("^([^ :]+)[ :]+([^ ]+)$", node_BS_sister)       ##################
                name_targetNode = match.group(1)
                bsBaclue = match.group(2)
                duplicationStatus = match.group(3)
                line_FN += "BS_of_" + name_targetNode + "_monophyly@" + bsBaclue + ","
                line_FN += "dupStatus_" + name_targetNode + "@" + duplicationStatus + ","
        else:
            #for targetSpeciesNode in speciesNodes_including_querySpecies:
            for targetSpeciesNode in childSpeciesNodes_orthogroup_including_querySpecies:
                name_targetSpeciesNode = make_nodeName_from_nodeLavel_NHXstyle(targetSpeciesNode[2])
                #print("name_targetSpeciesNode", name_targetSpeciesNode)
                line_FN += "BS_of_" + name_targetSpeciesNode + "_monophyly@NONE" + ","
                line_FN += "dupStatus_" + name_targetSpeciesNode + "@NONE" + ","

        if ">SisterGeneGroups" in  seqDict.keys():
            #print("seqDict[>Number_of_duplicatedNode]", seqDict[">Number_of_duplicatedNode"])
            for node_BS_sister in seqDict[">SisterGeneGroups"]:
                #print("node_BS_sister", node_BS_sister, "|")
                match = re.search("^([^ ]+) +([^ ]+) +([^ ]+)$", node_BS_sister)   ##################
                #match = re.search("^([^ :]+)[ :]+([^ ]+) +([^ ]+)$", node_BS_sister)   ##################
                name_targetNode = match.group(1)
                bsBaclue = match.group(2)
                name_sisterNode = match.group(3)
                line_FN += "Sister_of_" + name_targetNode + "@" + name_sisterNode + ","
                line_FN += "BS_with_" + name_targetNode + "@" + bsBaclue + ","
        else:
            #for targetSpeciesNode in speciesNodes_including_querySpecies:
            for targetSpeciesNode in childSpeciesNodes_orthogroup_including_querySpecies:
                name_targetSpeciesNode = make_nodeName_from_nodeLavel_NHXstyle(targetSpeciesNode[2])
                #print("name_targetSpeciesNode", name_targetSpeciesNode)
                line_FN += "Sister_of_" + name_targetSpeciesNode + "@NONE" + ","
                line_FN += "BS_with_" + name_targetSpeciesNode + "@NONE" + ","

        #if ">Number_of_duplicatedNode" in  seqDict.keys():
        #    #print("seqDict[>Number_of_duplicatedNode]", seqDict[">Number_of_duplicatedNode"])
        #    for node_dup in seqDict[">Number_of_duplicatedNode"]:
        #        #print("node_dup", node_dup)
        #        match = re.search("^([^ ]+) +(\d+)$", node_dup)
        #        node_name = match.group(1)
        #        num_dup = match.group(2)
        #        line_FN += "dup_" + node_name + "@" + num_dup + ","
        #else:
        #    #for targetSpeciesNode in speciesNodes_including_querySpecies:
        #    for childSpeciesNode_of_orthogroup in childSpeciesNodes_orthogorup:
        #        if not name_querySpecies in childSpeciesNode_of_orthogroup[1]:
        #            continue
        #        speciesNodeName = make_nodeName_from_nodeLavel_NHXstyle(childSpeciesNode_of_orthogroup[2])
        #        #print("name_targetSpeciesNode", name_targetSpeciesNode)
        #        line_FN += "dup_" + speciesNodeName + "@" + ","


        line_FN_lastComman_deleted = line_FN[:-1]
        lines_FN.append(line_FN_lastComman_deleted)

    return lines_FN


def print_csv(lines):
    out = open("results.csv", "w")
    flag = 0

    indexLine = make_indexLine(lines[0])
    out.write(indexLine + "\n")

    for line in lines:
        if AddintHeaderAfterAT == "D":
            #print("line1", line)
            line = re.sub(",[^@]+@", ",", line)   ########################### comment out for make sure
            line = re.sub("^[^@]+@", "", line)    ########################### comment out for make sure
            #print("line2", line)
        out.write(line + "\n")
    out.close()


def make_indexLine(first_line):
    indexLine_FN = re.sub("@[^,]+,", ",", first_line)
    indexLine_FN = re.sub("@[^@]+$", "", indexLine_FN)
    indexLine_FN = re.sub("@,", ",", indexLine_FN)
    indexLine_FN = re.sub("@$", "", indexLine_FN)
    return indexLine_FN

def copy_alignment_orthogroup():
    if os.path.exists(eachDirAddress + "180_aln_nucl_fas.txt"):
        #print("Present")
        recsTMP = readFasta_dict(eachDirAddress, "180_aln_nucl_fas.txt")
        #print("queryID", queryID)
        out = open(alignment_orthogroups + "/" + queryID + ".txt", "w")
        for nameLine, seq in recsTMP.items():
            out.write(nameLine + "\n")
            out.write(seq + "\n")
        out.close()



###############################################################



######################################################################################################################
#################################### Main program ####################################################################
######################################################################################################################

startTime = time.time()


#####

dbLinesTMP, taxonSamplingListTMP, SpeciesTreeTMP, blastEvalue, \
Number_of_hits_to_report_per_genome, Aligned_site_rate, dataset, BSthreshold, \
treeSearchMethod, num_rootSequences, keyNode, name_querySpecies, queryDatabase, \
dbAddress, outdir, alignment_orthogroups, mode, Switch_deleteIntermediateFiles, speciesWithGeneFunction\
= read_controlFile()




check_mode()

eachDirAddress = ""
if mode == "E" or mode == "D":
    check_toolsDirectory()

    print("\n\n############### " + queryID + " ################\n\n")
    #print("sys.version", sys.version)
    
    
    if not os.path.exists(dbAddress):
        print("Chack your database address. Stop.")
        print("Exit.")
    eachDirAddress = outdir + "/" + queryID + "/"
    #outdir = path_currentDirectory + "/" + outdir
    
    dirFileMake(queryDatabase, name_querySpecies, queryID)

else:
    eachDirAddress = "./"


#### Species Tree raddrrize (Top left)
#print("eachDirAddress",eachDirAddress)
make_treeFile("000_speciesTreeTMP.txt", SpeciesTreeTMP)
laderrizedTree = "tools/Rscript scripts/ladderizeTree.R " + eachDirAddress + "000_speciesTreeTMP.txt " + eachDirAddress + "000_speciesTree.txt"
#print("laderrizedTree: ", laderrizedTree)
subprocess.call(laderrizedTree, shell=True)
fs_SpeciesTree = open(eachDirAddress + "000_speciesTree.txt")
SpeciesTree = list(fs_SpeciesTree)[0]
fs_SpeciesTree.close()
SpeciesTree = SpeciesTree.rstrip("\n")


#### Reorder dbLines taxonSamplingList along with species tree
dbLines, taxonSamplingList = reorder_dbLines(dbLinesTMP, taxonSamplingListTMP, name_querySpecies, "000_speciesTree.txt")


allNodes_speciesTree  = collect_nodes_from_speciesTree()

if mode == "D":
    print("##### Tree draw ######")
    #print("queryID", queryID)
    if not os.path.exists(eachDirAddress + "100_analysisSummary.txt"):
        print("Error: cannot find 100_analysisSummary.txt in ", eachDirAddress)
        print("Your mode is:", mode)
        print("Mode D can be conducted after Mode E analyses.")
        print("So, choose Mode E if this is your first trial of ORTHOSCOPE*.")
        print("Also, make sure your command line.")
        print("Mode E or mode D allows only the geneID such as ENSORLT00000017526.1.")
        print("Your command line:")
        print(sys.argv[0], sys.argv[1])
        exit()

    print("##### 1st tree: APE (tree draw) ######")
    treePlotR_1st = "tools/Rscript scripts/treePlot.R " + eachDirAddress + "100_analysisSummary.txt " + " 1st_gene_tree_newick 1st_rearranged_gene_tree_newick Orthogroup " + eachDirAddress + "115_1st > " + eachDirAddress + "115_logTreePlotB.txt"
    #print("treePlotR: ", treePlotR_1st)
    subprocess.call(treePlotR_1st, shell=True)

    print("##### 2nd tree: APE (tree draw) ######")
    treePlotR = "tools/Rscript scripts/treePlot.R " + eachDirAddress + "100_analysisSummary.txt " + " 2nd_gene_tree_newick 2nd_rearranged_gene_tree_newick Rooting " + eachDirAddress + "240_2nd > " + eachDirAddress + "240_logTreePlotB.txt"
    #print("treePlotR: ", treePlotR)
    subprocess.call(treePlotR, shell=True)

    #print("outdir", outdir)
    #print("queryID", queryID)
    #print("    ",     )
    make_resHtml2(resHTMLlines_2steps)
    #make_resHtml_link_form_outside()
    exit()


if mode == "E":
    check_presense_of_databases()
    makeblastdb_database()

orthogroup_speciesNode = identifiy_orthogroup_speciesNode()

#querySpeciesNode = identify_speciesNode(name_querySpeciesNode)
querySpeciesNode = identify_speciesNode(name_querySpecies)
speciesNodes_including_querySpecies = collect_ancestralNodes(allNodes_speciesTree, querySpeciesNode)
#for targetSpeciesNode in speciesNodes_including_querySpecies:
#    name_targetSpeciesNode = make_nodeName_from_nodeLavel_NHXstyle(targetSpeciesNode[2])
#    print("name_targetSpeciesNode1", name_targetSpeciesNode)



childSpeciesNodes_orthogroup_including_querySpecies = collect_childNodesincluding_querySpecies(allNodes_speciesTree, orthogroup_speciesNode)
#for targetSpeciesNode in childSpeciesNodes_orthogroup_including_querySpecies:
#    name_targetSpeciesNode = make_nodeName_from_nodeLavel_NHXstyle(targetSpeciesNode[2])
#    print("name_targetSpeciesNode2", name_targetSpeciesNode)
#exit()


if mode == "S":
    #print("Mode S")
    outdir = outdir + "/"
    
    fileName_sum_list = sys.argv[1]
    geneIDsTMP = open(fileName_sum_list)
    geneIDs_TMP = list(geneIDsTMP)
    geneIDs = []
    for geneID in geneIDs_TMP:
        geneID = geneID.rstrip("\n")
        if not geneID:
            continue
        geneIDs.append(geneID)
    geneIDsTMP.close()
    
    lines_atmarkSeparated = make_lines_atmarkSeparated(geneIDs)
    #for line in lines_atmarkSeparated:
    #    print("line", line)
    #exit()
    print_csv(lines_atmarkSeparated)

    address_file = eachDirAddress + "000_speciesTreeTMP.txt"
    line_rm = "rm " + address_file
    #print("line_rm", line_rm)
    subprocess.call(line_rm, shell=True)
    exit()


if draw_speciesTree == "Draw":
    if not os.path.exists(outdir + "/speciesTree.pdf"):
        print("##### SpeciesTree draw ######")
        treePlot_speciesTree = "tools/Rscript scripts/treePlot.R control.txt " + " SpeciesTree dummy_rearranged_species_tree_newick Rooting speciesTree > speciesTree.pdf"
        #print("treePlot_speciesTree: ", treePlot_speciesTree)
        subprocess.call(treePlot_speciesTree, shell=True)
        new_path = shutil.move('./speciesTree.pdf', outdir)


checkUplodedFileAsFastaForamt()


'''


aaSeqMaker()

print("\n\n##### 1st tree: BLAST ######\n\n")
blastpSearch()


hitRecPicker()

print("\n\n##### 1st tree: MAFFT 1st round ######\n\n")
maffLine1 = "tools/mafft " + eachDirAddress + "030_retrievedAAfas.txt > " + eachDirAddress + "040_mafOutAA.txt"
print("maffLine1", maffLine1)
subprocess.call(maffLine1, shell=True)
#exit()
'''

fTMP = open(eachDirAddress + "040_mafOutAA.txt")
fMafOut = list(fTMP)
fTMP.close()
if not fMafOut:
    makeSummary(outfile_summary = "100_analysisSummary.txt")
    result = "No mafft out."
    print (result + " <br>")
    if Switch_deleteIntermediateFiles == "L":
        error_resHtmlMaker(result)
    if Switch_deleteIntermediateFiles == "D":
        deleteFiles()
    exit()


print("\n\n##### 1st tree: TRIMAL 1st round ######\n\n")
trimLine1 = "tools/trimal -out " + eachDirAddress + "042_AA.fas.trm -htmlout " + eachDirAddress + "042_AA.fas.trm.html -in " + eachDirAddress + "040_mafOutAA.txt -gappyout"
#print("trimLine1:", trimLine1, "\n");
subprocess.call(trimLine1, shell=True)

fTMP = open(eachDirAddress + "042_AA.fas.trm")
fTrimOut = list(fTMP)
fTMP.close()
if " 0 bp" in fTrimOut[0]:
     result = "Zero_trimledseq."
     error_makeSummary(result)
     print (result + " <br>")
     error_resHtmlMaker(result)
     exit()


print("\n\n##### 1st tree: calcilation for ShortSequence_threshold ######\n\n")
delete_sequences_with_alignedSiteRate()
#exit()


print("\n\n##### 1st tree: MAFFT 2nd round ######\n\n")
maffLine2 = "tools/mafft " + eachDirAddress + "044_overRateAA.fas > " + eachDirAddress + "050_retAAseqs.maf"
##print(maffLine2)
subprocess.call(maffLine2, shell=True)

print("\n\n##### 1st tree: TRIMAL 2nd round ######\n\n")
trimLine2 = "tools/trimal -out " + eachDirAddress + "052_AA.fas.trm -htmlout " + eachDirAddress + "052_AA.fas.trm.html -in " + eachDirAddress + "050_retAAseqs.maf -gappyout"
#print(trimLine2, "\n");
subprocess.call(trimLine2, shell=True)

print("\n\n##### 1st tree: PAL2NAL ######\n\n")
pal2nalLine = "tools/pal2nal.pl " + eachDirAddress + "050_retAAseqs.maf " + eachDirAddress + "044_overRateDNA.fas -output fasta > " + eachDirAddress +"054_p2nOutcDNAfas.txt"
#print (pal2nalLine)
subprocess.call(pal2nalLine, shell=True)
#print("<br>")
NJtreeFile = open(eachDirAddress + "054_p2nOutcDNAfas.txt")
NJtreeFileCont = list(NJtreeFile)
if not NJtreeFileCont:
    makeSummary(outfile_summary = "100_analysisSummary.txt")
    result = "No pal2nal out."
    print (result)
    #error_makeSummary(result)
    if Switch_deleteIntermediateFiles == "L":
        error_resHtmlMaker(result)
    else:
        deleteFiles()
    exit()

#trimaledv12_FileMakerDNA("054_p2nOutcDNAfas.txt", "052_AA.fas.trm.html", outfile="080_trimedCDNAPhy.txt")
trimaledv141_FileMakerDNA("054_p2nOutcDNAfas.txt", "052_AA.fas.trm.html", outfile="080_trimedCDNAPhy.txt")
#exit()

fas2phy(fastaFileName="052_AA.fas.trm", outPhyFileName="080_trimedAAPhy.txt")

print("\n\n##### 1st tree: APE (tree search) ######\n\n")
### NJ
outgroup1 = outGroupSelect("080_trimedAAPhy.txt")
if dataset == "Exclude3rd":
    print("The 1st gene tree is estimated by excluding 3rd codon positions.")
    #NJBSline1 = "tools/fastme -i " + eachDirAddress + "082_trimedBlockExc3rdFastmePhy.txt -d F84 -m BioNJ -b 100 -v 3 -o " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085_fastmelog.txt"
    phyCodonToBlock("080_trimedCDNAPhy.txt", 2, outfile="082_trimedBlockExc3rdPhy.txt")
    NJBSline1 = "tools/Rscript scripts/NJBS.R " + eachDirAddress + "082_trimedBlockExc3rdPhy.txt " + outgroup1 + " " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085_log.txt"
elif dataset == "Include3rd":
    #NJBSline1 = "scripts/fastme -i " + eachDirAddress + "082_trimedBlockInc3rdFastmePhy.txt -d F84 -m BioNJ -b 100 -v 3 -o " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085_fastmelog.txt"
    print("The 1st gene tree is estimated by including 3rd codon positions.")
    phyCodonToBlock("080_trimedCDNAPhy.txt", 3, outfile="082_trimedBlockInc3rdPhy.txt")
    NJBSline1 = "tools/Rscript scripts/NJBS.R " + eachDirAddress + "082_trimedBlockInc3rdPhy.txt " + outgroup1 + " " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085_log.txt"
else:
    print("The 1st gene tree is estimated using amino acid sequences.")
    phy2fastmePhy(phyFileName="080_trimedAAPhy.txt", outFastmePhyFileName="082_trimedAAFastmePhy.txt")
    NJBSline1 = "tools/fastme -i " + eachDirAddress + "082_trimedAAFastmePhy.txt --protein=WAG -m NJ -b 100 -v 3 -o " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085_fastmelog.txt"
#print("NJBSline1: ", NJBSline1)
subprocess.call(NJBSline1, shell=True)
#exit()

if not os.path.isfile(eachDirAddress + "085_NJBS1st.txt"):
    result = "No 1st tree. "
    #error_makeSummary(result)
    makeSummary(outfile_summary = "100_analysisSummary.txt")
    if Switch_deleteIntermediateFiles == "L":
        error_resHtmlMaker(result)
    if Switch_deleteIntermediateFiles == "D":
        deleteFiles()
    exit()

NJtreeFile = open(eachDirAddress + "085_NJBS1st.txt")
NJtreeFileCont = list(NJtreeFile)
if not NJtreeFileCont:
    result = "No 1st tree. "
    #error_makeSummary(result)
    makeSummary(outfile_summary = "100_analysisSummary.txt")
    if Switch_deleteIntermediateFiles == "L":
        error_resHtmlMaker(result)
    if Switch_deleteIntermediateFiles == "D":
        deleteFiles()
    exit()

print("\n\n##### 1st tree: NOTUNG ######\n\n")
NOTUNG1stLine = "java -jar tools/Notung.jar -s " + eachDirAddress + "000_speciesTree.txt -g " + eachDirAddress + "085_NJBS1st.txt --outputdir " + eachDirAddress + " --rearrange --threshold " + str(BSthreshold) + " --speciestag prefix  --maxtrees 5 --nolosses --treeoutput nhx > " + eachDirAddress + "085_NOTUNGlog.txt"
#print("NOTUNG1stLine:", NOTUNG1stLine)
subprocess.call(NOTUNG1stLine, shell=True)

if not os.path.isfile(eachDirAddress + "085_NJBS1st.txt.rearrange.0"):
    print ("Error in NOTUNG: Cannot compare the NJ and species tree.<br>")
    print ("Check your species tree.<br>")
    exit()



print("\n\n##### 1st tree: Making summary ######\n\n")
makeSummary(outfile_summary = "100_analysisSummary.txt")
#exit()


if Switch_deleteIntermediateFiles == "L":
    print("\n\n##### 1st tree: APE (tree draw) ######\n\n")
    treePlotR_1st = "tools/Rscript scripts/treePlot.R " + eachDirAddress + "100_analysisSummary.txt " + " 1st_gene_tree_newick 1st_rearranged_gene_tree_newick Orthogroup " + eachDirAddress + "115_1st > " + eachDirAddress + "115_logTreePlotB.txt"
    #print("treePlotR: ", treePlotR_1st)
    subprocess.call(treePlotR_1st, shell=True)
#### 1st tree search/rearrangemnet Finished
########################################################


########################################################
#### 2nd tree search/rearrangemnet Start
print("\n\n##### 2nd tree ######\n\n")


resDict_1st = readRes_dict(eachDirAddress + "100_analysisSummary.txt")
#print('resDict_1st[">BS_of_orthogroupBasalNode"][0]', resDict_1st[">BS_of_orthogroupBasalNode"][0])
#exit()
#if resDict_1st[">BS_of_orthogroupBasalNode"][0] == "No_orthogroup":
if resDict_1st[">BS_of_orthogroupBasalNode"][0].startswith("noOrthogroup_"):
    result = resDict_1st[">BS_of_orthogroupBasalNode"][0]
    #error_makeSummary(result)
    if Switch_deleteIntermediateFiles == "L":
        error_resHtmlMaker(result)
    if Switch_deleteIntermediateFiles == "D":
        deleteFiles()
    exit()

    #resHtmlMaker(resHTMLlines_2steps)
    if Switch_deleteIntermediateFiles == "D":
        deleteFiles()
    exit()
if len(resDict_1st[">Orthogroup"]) < 4:
    result = "Less than 4 orthogroup members."
    #error_makeSummary(result)
    #print (result + " <br>")
    if Switch_deleteIntermediateFiles == "L":
        error_resHtmlMaker(result)
    if Switch_deleteIntermediateFiles == "D":
        deleteFiles()
    exit()
    

rootLeaves = selectRootSp4secondTreeSearch()
make_2ndanalysis_seqFile(rootLeaves, outfile="150_noGapCDNAfas.txt")

cDNAfas2noGapAAFasFile("150_noGapCDNAfas.txt", outfile="150_noGapAA.txt")

print("\n\n##### 2nd tree: MAFFT ######\n\n")
maffLine2 = "mafft " + eachDirAddress + "/150_noGapAA.txt > " + eachDirAddress + "160_mafOut.txt"
##print("maffLine2", maffLine2)
subprocess.call(maffLine2, shell=True)

print("\n\n##### 2nd tree: TRIMAL ######\n\n")
trimLine2 = "tools/trimal -out " + eachDirAddress + "/170_trimedAAOutFas.txt -htmlout "  + eachDirAddress + "/170_aln_prot.html -in "  + eachDirAddress + "160_mafOut.txt -gappyout"
###print("trimLine2:", trimLine2, "\n");
subprocess.call(trimLine2, shell=True)

print("\n\n##### 2nd tree: PAL2NAL ######\n\n")
pal2nalLine = "tools/pal2nal.pl " + eachDirAddress + "160_mafOut.txt " + eachDirAddress + "150_noGapCDNAfas.txt -output fasta > " + eachDirAddress + "180_aln_nucl_fas.txt"
#print (pal2nalLine)
subprocess.call(pal2nalLine, shell=True)

fas2phy(fastaFileName = "160_mafOut.txt", outPhyFileName = "190_aln_prot.txt")
fas2phy(fastaFileName = "180_aln_nucl_fas.txt", outPhyFileName = "190_aln_nucl.txt")

#trimaledv12_FileMakerDNA("180_aln_nucl_fas.txt", "170_aln_prot.html", outfile = "200_trimedCDNAPhy.txt")
trimaledv141_FileMakerDNA("180_aln_nucl_fas.txt", "170_aln_prot.html", outfile = "200_trimedCDNAPhy.txt")
fas2phy("170_trimedAAOutFas.txt", "200_trimedAAPhy.txt")

phyCodonToBlock("200_trimedCDNAPhy.txt", 2, outfile="210_trimedBlockExc3rdPhy.txt")
phyCodonToBlock("200_trimedCDNAPhy.txt", 3, outfile="210_trimedBlockInc3rdPhy.txt")



#### 2nd tree search Start ####
print("\n\n##### 2nd tree: APE (tree search) ######\n\n")
outputRAxML = ""
#outGroup = outGroupSelect("210_trimedBlockExc3rdPhy.txt")
resDict_OG = readRes_dict(eachDirAddress + "100_analysisSummary.txt")
outGroup = resDict_OG[">Rooting"][0]
if treeSearchMethod == "ML":
    make_raxmlPartitionFile(outPartFile="220_raxmlPart.txt")
    raxmlLine = "scripts/raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -# 100 -m GTRCAT -s " + eachDirAddress + "210_trimedBlockExc3rdPhy.txt -w " + eachDirAddress + " -q " + eachDirAddress + "220_raxmlPart.txt -o " + outGroup + " -n txt -T 2"
    #print("raxmlLine:", raxmlLine)
    subprocess.call(raxmlLine, shell=True)
    moveRAxMLfiles(outfile = "230_2ndtree.txt")
    #exit()
else:
    if dataset == "Exclude3rd":
        print("The 2nd gene tree is estimated by excluding 3rd codon positions.")
        NJBSline2 = "tools/Rscript scripts/NJBS.R " + eachDirAddress + "210_trimedBlockExc3rdPhy.txt " + outGroup + " " + eachDirAddress + "230_2ndtree.txt > " + eachDirAddress + "230_log.txt"
    elif dataset == "Include3rd":
        print("The 2nd gene tree is estimated by including 3rd codon positions.")
        NJBSline2 = "tools/Rscript scripts/NJBS.R " + eachDirAddress + "210_trimedBlockInc3rdPhy.txt " + outGroup + " " + eachDirAddress + "230_2ndtree.txt > " + eachDirAddress + "230_log.txt"
    else:
        print("The 2nd gene tree is estimated using amino acid sequences. Stopped.")
        NJBSline2 = "tools/fastme -i " + eachDirAddress + "200_trimedAAPhy.txt --protein=WAG -m NJ -b 100 -v 3 -o " + eachDirAddress + "230_2ndtree.txt > " + eachDirAddress + "230_log.txt"
    #print("NJBSline2: ", NJBSline1)
    subprocess.call(NJBSline2, shell=True)


if not os.path.isfile(eachDirAddress + "230_2ndtree.txt"):
    print ("Error: Cannot estimate 2nd tree.<br>")
    #resHtmlMaker(resHTMLlines_2steps)
    if Switch_deleteIntermediateFiles == "D":
        deleteFiles()
    exit()

Tree2ndFile = open(eachDirAddress + "230_2ndtree.txt")
Tree2ndFileCont = list(Tree2ndFile)
if not Tree2ndFileCont:
    result = "2nd Tree Error: Cannot estimate the 2nd NJ tree."
    print (result + "<br>")
    make_resHtml2(resHTMLlines_2steps)
    if Switch_deleteIntermediateFiles == "D":
        deleteFiles()
    exit()
#### 2nd tree search End ####

rootBS100R = "tools/Rscript scripts/rootBS100.R " + eachDirAddress + "230_2ndtree.txt " + eachDirAddress + "230_2ndtreeRootBS100.txt"
#print("treePlotR: ", rootBS100R)
subprocess.call(rootBS100R, shell=True)

print("\n\n##### 2nd tree: NOTUNG ######\n\n")
NOTUNG2ndLine = "java -jar tools/Notung.jar -s " + eachDirAddress + "000_speciesTree.txt -g " + eachDirAddress + "230_2ndtreeRootBS100.txt --outputdir " + eachDirAddress + " --rearrange --threshold " + str(BSthreshold) + " --speciestag prefix  --maxtrees 5 --nolosses --treeoutput nhx > " + eachDirAddress + "230_NOTUNGlog.txt"
#print("NOTUNG2ndLine:", NOTUNG2ndLine)
subprocess.call(NOTUNG2ndLine, shell=True)
#exit()


print("\n\n##### 2nd tree: Making summary ######\n\n")
add_makeSummary(outfile_summary2 = "100_analysisSummary.txt")
#exit()


if Switch_deleteIntermediateFiles == "L":
    print("\n\n##### 2nd tree: APE (tree draw) ######\n\n")
    treePlotR = "tools/Rscript scripts/treePlot.R " + eachDirAddress + "100_analysisSummary.txt " + " 2nd_gene_tree_newick 2nd_rearranged_gene_tree_newick Rooting " + eachDirAddress + "240_2nd > " + eachDirAddress + "240_logTreePlotB.txt"
    #print("treePlotR: ", treePlotR)
    subprocess.call(treePlotR, shell=True)
    make_resHtml2(resHTMLlines_2steps)

copy_alignment_orthogroup()

if Switch_deleteIntermediateFiles == "D":
    deleteFiles()

exit()


