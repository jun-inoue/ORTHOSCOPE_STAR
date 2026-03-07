#!/usr/local/bin/python3

import sys
import re, os, shutil
from collections import OrderedDict
import subprocess
import time
import argparse

AddintHeaderAfterAT = "D"   ## L:leave or D:Delete @xxxx for the summarize analysis.
draw_speciesTree = "Not"  ## Draw: Draw species tree in the .pdf file. or Not:


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
        <td align="center" valign="top">NJ tree (<a href="EACHDIRADDRESS_240_2ndGeneTree.pdf" target="_blank">PDF</a>)</td>
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
class OrthoScopeContext:
    """解析全体で共有する設定値・中間生成物の置き場所（状態の箱）"""
    def __init__(self):
        # --- control.txt から取得される設定（read_controlFile の戻り） ---
        self.mode = None
        self.keyNode = None
        self.name_querySpecies = None
        self.dbAddress = None
        self.toolAddress = None
        self.scriptAddress = None
        self.outdir = None
        self.outdir_4_ReE1st = None
        self.alignment_orthogroups = None
        self.dataset = None
        self.BSthreshold = None
        self.BSthreshold_4_ReE1st = None
        self.treeSearchMethod = None
        self.num_rootSequences = None
        self.Switch_deleteIntermediateFiles = None
        self.speciesWithGeneFunction = None
        self.dbLines = None
        self.taxonSamplingList = None
        self.queryDatabase = None
        self.blastEvalue = None
        self.Number_of_hits_to_report_per_genome = None
        self.aligned_site_rate = None
        self.startTime = None

        # read_controlFile 後に reorder 前の一時保持
        self._dbLinesTMP = None
        self._taxonSamplingListTMP = None
        self._SpeciesTreeTMP = None

        # --- 解析の途中で生まれる中間生成物 ---
        self.queryID = None
        self.eachDirAddress = None
        self.eachDirAddress_e1stre = None
        self.SpeciesTree = None
        self.allNodes_speciesTree = None
        self.focalNode_speciesTree = None
        self.childSpeciesNodes_AllGroup = None
        self.childSpeciesNodes_focalGroup = None
        self.SpeciesTree = None
        self.dbLines = None
        self.taxonSamplingList = None

        self.allNodes_speciesTree = None
        self.focalNode_speciesTree = None
        self.speciesNodes_including_querySpecies = None
        self.childSpeciesNodes_AllGroup = None
        self.childSpeciesNodes_focalGroup = None

        self.resDict_1st = None

        # 他、必要に応じて追加していく（後で徐々に ctx に集約）


def run_mode_E_pipeline(ctx):

    # --- ReE1st: 既存1st-tree/SpeciesTreeでNotung再配置のみ（新規ディレクトリへ出力） ---

    if ctx.mode == "ReE1st":
        print("##### Mode ReE1st analysis ######")

        #print("ctx.eachDirAddress", ctx.eachDirAddress)
        #print("### line 448")
        #exit()

        os.makedirs(ctx.eachDirAddress_e1stre, exist_ok=True)

        _e_tree_rearrange_e1stre(ctx)

        # （任意）図の描画とHTML
        if ctx.Switch_deleteIntermediateFiles == "L":
            if ctx.BSthreshold_4_ReE1st == "reconcile":
                tasked_tree_newick = "1st_reconciled_gene_tree_newick"
            else:
                tasked_tree_newick = "1st_rearranged_gene_tree_newick"
            treePlotR_2ReE1st = (
                f"{ctx.toolAddress}Rscript {ctx.scriptAddress}treePlot.R "
                f"{ctx.eachDirAddress_e1stre}100_analysisSummary.txt "
                f"1st_gene_tree_newick "
                f"{tasked_tree_newick} "
                f"Rooting_4_1stTree Rooting_4_1stTree "
                f"{ctx.eachDirAddress_e1stre}115_1st > {ctx.eachDirAddress_e1stre}115_logTreePlotB.txt"
            )
            print("treePlotR_2ReE1st:", treePlotR_2ReE1st)
            subprocess.call(treePlotR_2ReE1st, shell=True)

            make_resHtml2(
                ctx.queryID,
                ctx.mode,
                ctx.BSthreshold_4_ReE1st,
                ctx.eachDirAddress_e1stre,
                resHTMLlines_2steps
            )

        if ctx.Switch_deleteIntermediateFiles == "D":
            deleteFiles(ctx.eachDirAddress_e1stre)

        elapsed = time.time() - ctx.startTime
        print(f"#####")
        print(f"Mode {ctx.mode} analysis has completed successfully. Elapsed time: {elapsed:.1f} seconds.")
        print("No 2nd-tree estimation is performed for mode ReE1st.\n")
        exit()

    # --- 既存の E/E1st 経路 ---

    _e_prepare_databases(ctx)

    _e_run_1st_tree_blast_to_trim(ctx)
    _e_run_1st_tree_search_notung_summary(ctx)

    # E1st/ReE1st はここで終了
    if ctx.mode in ("E1st", "ReE1st"):
        if ctx.Switch_deleteIntermediateFiles == "D":
            deleteFiles(ctx.eachDirAddress)
        elapsed = time.time() - ctx.startTime
        print(f"#####")
        print(f"Mode {ctx.mode} analysis has completed successfully. Elapsed time: {elapsed:.1f} seconds.")
        print("No 2nd-tree estimation is performed for mode E1st or ReE1st. \n\n")
        exit()

    _e_run_2nd_tree_pipeline(ctx)
    _e_finalize(ctx)

def run_or_make_summary_and_exit(ctx, cmd: str, reason: str, exit_code: int = 0):
    """
    外部コマンドを実行し、失敗したら 100_analysisSummary.txt を必ず作って終了する。
    reason は "mafft (1st round) failed" のように短く書き、Summary には
    "FATAL_ERROR:<reason>" 形式で残す。
    """
    rc = subprocess.call(cmd, shell=True)
    if rc == 0:
        return 0

    fatal = f"FATAL_ERROR:{reason}"
    print(f"[ERROR] External command failed (exit={rc}):\n{cmd}")
    print(f"[ERROR] {fatal}")

    makeSummary(
        ctx.SpeciesTree, ctx.taxonSamplingList, ctx.mode, ctx.dataset,
        ctx.keyNode, ctx.startTime, ctx.eachDirAddress, ctx.eachDirAddress_e1stre, ctx.BSthreshold, ctx.BSthreshold_4_ReE1st,
        ctx.num_rootSequences, ctx.allNodes_speciesTree,
        ctx.childSpeciesNodes_AllGroup,
        aligned_site_rate=ctx.aligned_site_rate,
        outfile_summary="100_analysisSummary.txt",
        fatal_error_msg=fatal,
    )

    if ctx.Switch_deleteIntermediateFiles == "D":
        deleteFiles(ctx.eachDirAddress)

    sys.exit(exit_code)

#def run_or_die(cmd: str):
#    rc = subprocess.call(cmd, shell=True)
#    if rc != 0:
#        print(f"[ERROR] External command failed (exit={rc}):\n{cmd}")
#        sys.exit(rc)
#    return rc


def log(*args):
    # 将来 logging に差し替えやすい軽いラッパ
    print(*args)

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="orthoscope_star_v1.2.5.py",
        description="ORTHOSCOPE* main entry point",
        add_help=True,
    )
    parser.add_argument(
        "query_id",
        help="Mode E/E1st/ReE1st/D/D1st/reD1st は gene ID、Mode S/S1st/ReS1st は .txt（IDのリスト）"
    )
    return parser.parse_args(argv)

def initialize_context(query_id):
    # control.txt を読み、ctx を初期化
    (
        dbLinesTMP, 
        taxonSamplingListTMP, 
        SpeciesTreeTMP,
        blastEvalue, 
        Number_of_hits_to_report_per_genome, 
        aligned_site_rate,
        dataset, 
        BSthreshold, 
        BSthreshold_4_ReE1st, 
        treeSearchMethod, 
        num_rootSequences,
        keyNode, 
        name_querySpecies, 
        queryDatabase, 
        dbAddress, 
        toolAddress,
        scriptAddress, 
        outdir, 
        outdir_4_ReE1st, 
        alignment_orthogroups, 
        mode,
        Switch_deleteIntermediateFiles, 
        speciesWithGeneFunction
    ) = read_controlFile()
    
    ctx = OrthoScopeContext()
    # control.txt 由来
    ctx._dbLinesTMP = dbLinesTMP
    ctx._taxonSamplingListTMP = taxonSamplingListTMP
    ctx._SpeciesTreeTMP = SpeciesTreeTMP
    ctx.blastEvalue = blastEvalue
    ctx.Number_of_hits_to_report_per_genome = Number_of_hits_to_report_per_genome
    ctx.aligned_site_rate = aligned_site_rate
    ctx.dataset = dataset
    
    try:
        ctx.BSthreshold = int(BSthreshold)
    except ValueError:
        ctx.BSthreshold = BSthreshold.strip()
    try:
        ctx.BSthreshold_4_ReE1st = int(BSthreshold_4_ReE1st)
    except ValueError:
        ctx.BSthreshold_4_ReE1st = BSthreshold_4_ReE1st.strip()

    ctx.treeSearchMethod = treeSearchMethod
    ctx.num_rootSequences = num_rootSequences
    ctx.keyNode = keyNode
    ctx.name_querySpecies = name_querySpecies
    ctx.queryDatabase = queryDatabase
    ctx.dbAddress = dbAddress
    ctx.toolAddress = toolAddress
    ctx.scriptAddress = scriptAddress
    ctx.outdir = outdir
    ctx.outdir_4_ReE1st = outdir_4_ReE1st
    ctx.alignment_orthogroups = alignment_orthogroups
    ctx.mode = mode
    ctx.Switch_deleteIntermediateFiles = Switch_deleteIntermediateFiles
    ctx.speciesWithGeneFunction = speciesWithGeneFunction
    
    # 実行引数
    ctx.queryID = query_id
    ctx.startTime = time.time()

    if ctx.mode in ("E", "E1st", "D", "D1st"):
        #ctx.eachDirAddress = ctx.outdir + "/" + ctx.queryID + "/"
        ctx.eachDirAddress = os.path.join(ctx.outdir, ctx.queryID, "")
    elif ctx.mode == "ReE1st":
        ctx.eachDirAddress = os.path.join(ctx.outdir, ctx.queryID, "")
        ctx.eachDirAddress_e1stre = os.path.join(ctx.outdir_4_ReE1st, ctx.queryID, "")
    elif ctx.mode == "ReD1st":
        ctx.eachDirAddress_e1stre = os.path.join(ctx.outdir_4_ReE1st, ctx.queryID, "")

    return ctx

def echo_major_settings(ctx):
    control_path = os.path.abspath("control.txt")
    log("#### Major settings")
    log(f"{control_path} is used.")
    if not os.path.exists(control_path):
        log("WARNING: control.txt が見つかりません。スクリプト配置/実行ディレクトリを確認してください。")
    for k, v in (
        (">Mode", ctx.mode),
        (">QuerySpecies", ctx.name_querySpecies),
        (">Database", ctx.dbAddress),
        (">Outdir", ctx.outdir),
        (">Outdir_4_ReE1st", ctx.outdir_4_ReE1st),
        (">Dataset", ctx.dataset),
        (">BSthreshold", ctx.BSthreshold),
        (">BSthreshold_4_ReE1st", ctx.BSthreshold_4_ReE1st),
        (">Switch_deleteIntermediateFiles", ctx.Switch_deleteIntermediateFiles),
    ):
        log(k); log(v)
    if ctx.mode in ("E", "E1st", "ReE1st", "D", "D1st", "ReD1st"):
        log(f"############### {ctx.queryID} ################")


def prepare_files_and_species_tree(ctx):
    #print("### prepare_files_and_species_tree() ###")
    #print("ctx.eachDirAddress", ctx.eachDirAddress)
    #print("ctx.eachDirAddress_e1stre", ctx.eachDirAddress_e1stre)
    #exit()
    
    check_mode(ctx.mode, ctx.queryID)
    check_toolsDirectory(ctx.mode, ctx.toolAddress, ctx.dataset)
    check_scriptsDirectory(ctx.scriptAddress)
    
    #print("ctx.eachDirAddress", ctx.eachDirAddress, "|")
    #exit()

    # ---- D / D1st は描画専用モード：種系統樹ベースの前処理はスキップ ----
    if ctx.mode in ("D", "D1st"):
        summary = os.path.join(ctx.eachDirAddress, "100_analysisSummary.txt")
        if (not os.path.isfile(summary)) or (os.path.getsize(summary) == 0):
            log("Error: cannot find 100_analysisSummary.txt in ", ctx.eachDirAddress)
            sys.exit(1)
        return
    if ctx.mode == "ReD1st":
        summary = os.path.join(ctx.eachDirAddress_e1stre, "100_analysisSummary.txt")
        if (not os.path.isfile(summary)) or (os.path.getsize(summary) == 0):
            log("Error: cannot find 100_analysisSummary.txt in ", ctx.eachDirAddress_e1stre)
            sys.exit(1)
        return
    if ctx.mode in ("E", "E1st", "ReE1st"):
        if not os.path.exists(ctx.dbAddress):
            log("Chack your database address. Stop.")
            sys.exit(1)

    if ctx.mode in ("E", "E1st", "ReE1st"):
        dirFileMake(ctx.mode, ctx.outdir, ctx.eachDirAddress, ctx.alignment_orthogroups)

    if ctx.mode in ("E", "E1st"):
        make_querySeqFile(ctx.queryDatabase, ctx.name_querySpecies, ctx.dbAddress, ctx.eachDirAddress, ctx.queryID)
        check_uploaded_file_as_fasta_format(ctx.eachDirAddress)

    # species tree トップレフト
    in_tree = "./000_speciesTreeTMP.txt"
    out_tree = "./000_speciesTree_topLeft.txt"
    with open(in_tree, "w") as fhtree:
        fhtree.write(ctx._SpeciesTreeTMP + "\n")
    tree_topLeft_cmd = (
        f"{ctx.toolAddress}Rscript {ctx.scriptAddress}ladderizeTree.R "
        f"{in_tree} "
        f"{out_tree}"
    )
    subprocess.call(tree_topLeft_cmd, shell=True)
    if (not os.path.isfile(out_tree)) or (os.path.getsize(out_tree) == 0):
        print("[ERROR] ladderizeTree.R failed to create:", out_tree)
        sys.exit(1)

    with open("000_speciesTree_topLeft.txt") as f:
        ctx.SpeciesTree = f.readline().rstrip("\n")

    ctx.dbLines, ctx.taxonSamplingList = reorder_dbLines(
        ctx._dbLinesTMP, ctx._taxonSamplingListTMP, ctx.name_querySpecies,
        "./", "000_speciesTree_topLeft.txt"
    )
    

    (
        ctx.allNodes_speciesTree,
        ctx.focalNode_speciesTree,
        ctx.speciesNodes_including_querySpecies,
        ctx.childSpeciesNodes_AllGroup,
        ctx.childSpeciesNodes_focalGroup,
    ) = get_species_tree_info(
        ctx.SpeciesTree, ctx.keyNode, ctx.name_querySpecies
    )

def run_mode_S_block(ctx):
    log("### Mode S/S1st/ReS1st ###")
    if ctx.mode in ("S", "S1st"):
        outdir_S = ctx.outdir
    elif ctx.mode == "ReS1st":
        outdir_S = ctx.outdir_4_ReE1st
    else:
        print("Error in run_mode_S_block(). Check your mode ", ctx.mode)
    fileName_sum_list = ctx.queryID

    # --- ここを書き換え（for文版） ---
    geneIDs = []
    with open(fileName_sum_list) as f:
        for x in f:
            line = x.strip()
            if line:  # 空行を除外
                geneIDs.append(line)
    # --- ここまで ---

    lines_atmarkSeparated = make_lines_atmarkSeparated(
        mode=ctx.mode,
        SpeciesTree=ctx.SpeciesTree,
        keyNode=ctx.keyNode,
        allNodes_speciesTree=ctx.allNodes_speciesTree,
        childSpeciesNodes_AllGroup=ctx.childSpeciesNodes_AllGroup,
        childSpeciesNodes_focalGroup=ctx.childSpeciesNodes_focalGroup,
        speciesWithGeneFunction=ctx.speciesWithGeneFunction,
        outdir=outdir_S,
        geneIDs_fn=geneIDs,
    )
    print_csv(lines_atmarkSeparated)

    #address_file = "./" + "000_speciesTreeTMP.txt"
    #subprocess.call("rm " + address_file, shell=True)
    sys.exit(0)

def run_mode_D_draw_only(ctx):
    #print("### run_mode_D_draw_only() ###")
    #print("ctx.mode", ctx.mode)
    #exit()
    log("##### Tree draw ######")
    if ctx.mode in ("D", "D1st"):
        if not os.path.exists(ctx.eachDirAddress + "100_analysisSummary.txt"):
            log("Error: cannot find 100_analysisSummary.txt in ", ctx.eachDirAddress)
            log("Your mode is:", ctx.mode)
            log(sys.argv[0], ctx.queryID)
            sys.exit(1)
    elif ctx.mode == "ReD1st":
        if not os.path.exists(ctx.eachDirAddress_e1stre + "100_analysisSummary.txt"):
            log("Error: cannot find 100_analysisSummary.txt in ", ctx.eachDirAddress_e1stre)
            log("Your mode is:", ctx.mode)
            log(sys.argv[0], ctx.queryID)
            sys.exit(1)
    else:
        print("Error in run_mode_D_draw_only(): check your mode:", ctx.mode)
        exit()

    if ctx.mode in ("D", "D1st"):
        resDict_1stSummary = readRes_dict(ctx.eachDirAddress + "100_analysisSummary.txt")
        BS_of_orthogroupBasalNode = resDict_1stSummary[">BS_of_orthogroupBasalNode"][0]

    print("##### 1st tree: APE (tree draw) ######")
    if ctx.mode == "D1st":
        treePlotR_D1st = (
            f"{ctx.toolAddress}Rscript {ctx.scriptAddress}treePlot.R "
            f"{ctx.eachDirAddress}100_analysisSummary.txt "   #name_summaryFile
            f"1st_gene_tree_newick "                          #Gene_tree_newick
            f"1st_rearranged_gene_tree_newick "               #Rearranged_gene_tree_newick
            f"Rooting_4_1stTree "                             #rthick_branch_species_name_line
            f"Rooting_4_1stTree "                             #root_species_name_line Rooting_4_2ndTree Rooting_4_1stTree Orthogroup
            f"{ctx.eachDirAddress}115_1st > {ctx.eachDirAddress}115_logTreePlotB.txt"
        )
        print("treePlotR_D1st: ", treePlotR_D1st)
        subprocess.call(treePlotR_D1st, shell=True)
    elif ctx.mode == "ReD1st":
        if ctx.BSthreshold_4_ReE1st == "reconcile":
            tasked_tree_newick = "1st_reconciled_gene_tree_newick"
        else:
            tasked_tree_newick = "1st_rearranged_gene_tree_newick"
        treePlotR_D1st = (
            f"{ctx.toolAddress}Rscript {ctx.scriptAddress}treePlot.R "
            f"{ctx.eachDirAddress_e1stre}100_analysisSummary.txt "   #name_summaryFile
            f"1st_gene_tree_newick "                          #Gene_tree_newick
            f"{tasked_tree_newick} "               #Rearranged_gene_tree_newick
            f"Rooting_4_1stTree "                             #rthick_branch_species_name_line
            f"Rooting_4_1stTree "                             #root_species_name_line Rooting_4_2ndTree Rooting_4_1stTree Orthogroup
            f"{ctx.eachDirAddress_e1stre}115_1st > {ctx.eachDirAddress_e1stre}115_logTreePlotB.txt"
        )
        print("treePlotR_D1st: ", treePlotR_D1st)
        subprocess.call(treePlotR_D1st, shell=True)
    elif ctx.mode == "D":
        treePlotR_1st_D = (
            f"{ctx.toolAddress}Rscript {ctx.scriptAddress}treePlot.R "        
            f"{ctx.eachDirAddress}100_analysisSummary.txt "   #name_summaryFile
            f"1st_gene_tree_newick "                          #Gene_tree_newick
            f"1st_rearranged_gene_tree_newick "               #Rearranged_gene_tree_newick
            f"Orthogroup "                                    #thick_branch_species_name_line
            f"Rooting_4_2ndTree "                             #root_species_name_line Rooting_4_2ndTree Rooting_4_1stTree Orthogroup
            f"{ctx.eachDirAddress}115_1st > {ctx.eachDirAddress}115_logTreePlotB.txt"
        )
        print("treePlotR_1st_D: ", treePlotR_1st_D)
        subprocess.call(treePlotR_1st_D, shell=True)
    
        print("##### 2nd tree: APE (tree draw) ######")
        treePlotR_2nd_D = (
            f"{ctx.toolAddress}Rscript {ctx.scriptAddress}treePlot.R "
            f"{ctx.eachDirAddress}100_analysisSummary.txt " #name_summaryFile
            f"2nd_gene_tree_newick "                        #Gene_tree_newick
            f"2nd_rearranged_gene_tree_newick "             #Rearranged_gene_tree_newick
            f"Rooting_4_2ndTree "                           #thick_branch_species_name_line
            f"Rooting_4_2ndTree "                           #root_species_name_line Rooting_4_2ndTree Rooting_4_1stTree Orthogroup
            f"{ctx.eachDirAddress}240_2nd > {ctx.eachDirAddress}240_logTreePlotB.txt"
        )
        print("treePlotR_2nd_D: ", treePlotR_2nd_D)
        subprocess.call(treePlotR_2nd_D, shell=True)
    else:
        print("Error. mode shoud be D, D1st, ReD1st. Your mode is ", ctx.mode)
        exit()

    #print("outdir", outdir)
    #print("queryID", queryID)
    #print("    ",     )

    if ctx.mode in ("D", "D1st"):
        make_resHtml2(
            ctx.queryID,
            ctx.mode,
            ctx.BSthreshold_4_ReE1st,
            ctx.eachDirAddress,
            resHTMLlines_2steps
        )
    elif ctx.mode == "ReD1st":
        make_resHtml2(
            ctx.queryID,
            ctx.mode,
            ctx.BSthreshold_4_ReE1st,
            ctx.eachDirAddress_e1stre,
            resHTMLlines_2steps
        )
    else:
        print("Error2. mode shoud be D, D1st, ReD1st. Your mode is ", ctx.mode)
        exit()

    sys.exit(0)

def _e_prepare_databases(ctx):
    # 1) 事前チェック
    check_presense_of_databases(ctx.dbAddress, ctx.dbLines)
    makeblastdb_database(ctx.dbAddress, ctx.dbLines, ctx.toolAddress)

    # 2) 1st tree 前半の準備
    aaSeqMaker(ctx.eachDirAddress)
    

def _e_tree_rearrange_e1stre(ctx):
    #print("### _e_tree_rearrange_e1stre() ###")
    #print("ctx.eachDirAddress", ctx.eachDirAddress)
    #exit()
    """
    ReE1st:
      - ctx.eachDirAddress/queryID/100_analysisSummary.txt（既存）を読み、
        >1st_gene_tree_newick と >SpeciesTree を取得
      - その組で Notung を再実行（再配置）
      - 生成された木で
          >1st_rearranged_gene_tree_newick
          >1st_rearranged_gene_tree_NHX
        のセクションだけ差し替えて保存（保存先は ctx.eachDirAddress_e1stre 側）
      - 2nd-tree は実施しない
    """

    # 1) 入力は「元」の eachDirAddress 側
    src_summary = os.path.join(ctx.eachDirAddress, "100_analysisSummary.txt")
    if (not os.path.isfile(src_summary)) or (os.path.getsize(src_summary) == 0):
        print("[ERROR] ReE1st requires existing 100_analysisSummary.txt in:", ctx.eachDirAddress)
        print("       Run mode E/E1st beforehand, or place the file manually.")
        sys.exit(1)

    resDict_src = readRes_dict(src_summary)

    resDict_1stSummary = readRes_dict(ctx.eachDirAddress + "100_analysisSummary.txt")
    topHitName_1stQuery = resDict_1stSummary[">QuerySequence"][0]
    topHitName_1stQuery = re.sub(r" +.*", "", topHitName_1stQuery)


    # 2) 以降の生成物は「新しい」ctx.eachDirAddress_e1stre 側へ
    if ">1st_gene_tree_newick" not in resDict_src or not resDict_src[">1st_gene_tree_newick"]:
        print("[ERROR] >1st_gene_tree_newick not found in existing summary.")
        sys.exit(1)

    species_tree_newick = resDict_src[">SpeciesTree"][0]
    file_species_tree = os.path.join(ctx.eachDirAddress_e1stre, "000_speciesTree.txt")
    with open(file_species_tree, "w") as fhtree:
        fhtree.write(species_tree_newick.strip() + "\n")

    gene_tree_newick = resDict_src[">1st_gene_tree_newick"][0]
    file_gene_tree_newick = os.path.join(ctx.eachDirAddress_e1stre, "085_NJBS1st.txt")
    with open(file_gene_tree_newick, "w") as fhtree:
        fhtree.write(gene_tree_newick.strip() + "\n")

    # 3) Notung を再実行（再配置）
    print("##### 1st tree: NOTUNG (ReE1st) ######")
    #print("ctx.BSthreshold_4_ReE1st", ctx.BSthreshold_4_ReE1st)
    # BSthreshold_4_ReE1st が 0〜100 の数値かどうか判定
    if isinstance(ctx.BSthreshold_4_ReE1st, (int, float)) and 0 <= ctx.BSthreshold_4_ReE1st <= 100:
        notung_cmd = (
            f"java -jar {ctx.toolAddress}Notung.jar "
            f"-s {file_species_tree} "
            f"-g {file_gene_tree_newick} --outputdir {ctx.eachDirAddress_e1stre} "
            f"--rearrange --threshold {ctx.BSthreshold_4_ReE1st} --speciestag prefix "
            f"--maxtrees 5 --nolosses --treeoutput nhx > {ctx.eachDirAddress_e1stre}085_NOTUNGlog.txt"
        )
    elif ctx.BSthreshold_4_ReE1st == "reconcile":
        notung_cmd = (
            f"java -jar {ctx.toolAddress}Notung.jar "
            f"-s {file_species_tree} "
            f"-g {file_gene_tree_newick} --outputdir {ctx.eachDirAddress_e1stre} "
            f"--reconcile --speciestag prefix "
            f"--maxtrees 5 --nolosses --treeoutput nhx "
            f"> {ctx.eachDirAddress_e1stre}085_NOTUNGlog.txt"
        )
    elif ctx.BSthreshold_4_ReE1st == "transfer":
        notung_cmd = (
            f"java -jar {ctx.toolAddress}Notung.jar "
            f"-s {file_species_tree} "
            f"-g {file_gene_tree_newick} --outputdir {ctx.eachDirAddress_e1stre} "
            f"--reconcile --speciestag prefix --infertransfers true "
            f"--multsols 5 --nolosses --treeoutput nhx "
            f"> {ctx.eachDirAddress_e1stre}085_NOTUNGlog.txt"
        )
    else:
        print("Error. check task.")
        exit()
    #print("NOTUNG1stLine (ReE1st):", notung_cmd)
    #print("BSthreshold_4_ReE1st:", ctx.BSthreshold_4_ReE1st)
    #print("")
    subprocess.call(notung_cmd, shell=True)
    
    # ファイル名：reconcile を rearrange に書き換え
    if ctx.BSthreshold_4_ReE1st == "reconcile":
        src = f"{ctx.eachDirAddress_e1stre}085_NJBS1st.txt.reconciled"
        dst = f"{ctx.eachDirAddress_e1stre}085_NJBS1st.txt.rearrange.0"
        os.rename(src, dst)

    # rearrange/reconcile 結果を、それぞれ変数に格納
    file_rearranged_gene_tree_newick  = os.path.join(ctx.eachDirAddress_e1stre, "085_NJBS1st.txt.rearrange.0")
    if not os.path.isfile(file_rearranged_gene_tree_newick):
        print("[ERROR] Notung did not produce:", file_rearranged_gene_tree_newick)
        print("        Check your species tree / 1st gene tree.")
        sys.exit(1)
    with open(file_rearranged_gene_tree_newick) as f:
        rearranged_1st_gene_tree_NHX = f.readline().rstrip("\n")
    allGeneNodesSR_1stTree = collect_nodes_from_NHX(ctx.keyNode, rearranged_1st_gene_tree_NHX)

    # NHX → newick へ変換
    rearranged_tree_newick = change_nhx_to_newick_with_NHXnodeName(rearranged_1st_gene_tree_NHX)


    # 4) Summary で、該当 セクションだけ差し替え（他は保持）
    dst_summary = os.path.join(ctx.eachDirAddress_e1stre, "100_analysisSummary.txt")
    with open(dst_summary, "w") as out:
        for name, list_content in resDict_src.items():
            if name == ">1st_rearranged_gene_tree_newick":
                if ctx.BSthreshold_4_ReE1st == "reconcile":
                    out.write(">1st_reconciled_gene_tree_newick\n")
                else:
                   out.write(name + "\n")
                out.write(rearranged_tree_newick + "\n\n")
            elif name == ">1st_rearranged_gene_tree_NHX":
                if ctx.BSthreshold_4_ReE1st == "reconcile":
                    out.write(">1st_reconciled_gene_tree_NHX\n")
                else:
                    out.write(name + "\n")
                out.write(rearranged_1st_gene_tree_NHX + "\n\n")
            elif name == ">1st_gene_tree_newick":
                out.write(name + "\n")
                for line in list_content:
                    out.write(line + "\n")
                out.write("\n")

                out.write(">MonophyleticGeneGroups_1stTree\n")
                list_resLines_mono = make_list_resLines_monophyletic(ctx.allNodes_speciesTree, allGeneNodesSR_1stTree, ctx.childSpeciesNodes_AllGroup, topHitName_1stQuery)
                for line in list_resLines_mono:
                    out.write(line + "\n")
                out.write("\n")

                out.write(">SisterGeneGroups_1stTree\n")
                list_resLines_sister = make_list_resLines_sister(ctx.allNodes_speciesTree, allGeneNodesSR_1stTree, ctx.childSpeciesNodes_AllGroup, topHitName_1stQuery)
                for line in list_resLines_sister:
                    out.write(line + "\n")
                out.write("\n")

            elif name in (">MonophyleticGeneGroups_1stTree", ">SisterGeneGroups_1stTree"):
                pass
            else:
                out.write(name + "\n")
                for line in list_content:
                    out.write(line + "\n")
                out.write("\n")


def _e_run_1st_tree_blast_to_trim(ctx):
    print("\n\n##### 1st tree: BLAST ######\n\n")
    blastpSearch(
        ctx.dbLines, ctx.dbAddress, ctx.toolAddress,
        ctx.blastEvalue, ctx.Number_of_hits_to_report_per_genome,
        ctx.eachDirAddress
    )

    hitRecPicker(ctx.dbAddress, ctx.dbLines, ctx.eachDirAddress)

    print("\n\n##### 1st tree: MAFFT 1st round ######\n\n")
    maffLine1 = (
        f"{ctx.toolAddress}mafft "
        f"{ctx.eachDirAddress}030_retrievedAAfas.txt "
        f"> {ctx.eachDirAddress}040_mafOutAA.txt"
    )
    #subprocess.call(maffLine1, shell=True)
    run_or_make_summary_and_exit(ctx, maffLine1, "ERROR: mafft (1st round) failed")

    fTMP = open(ctx.eachDirAddress + "040_mafOutAA.txt")
    fMafOut = list(fTMP)
    fTMP.close()
    if not fMafOut:
        makeSummary(
            ctx.SpeciesTree, ctx.taxonSamplingList, ctx.mode, ctx.dataset,
            ctx.keyNode, ctx.startTime, ctx.eachDirAddress, ctx.eachDirAddress_e1stre, ctx.BSthreshold,
            ctx.num_rootSequences, ctx.allNodes_speciesTree,
            ctx.childSpeciesNodes_AllGroup,
            aligned_site_rate=ctx.aligned_site_rate,
            outfile_summary="100_analysisSummary.txt",
        )
        result = "No mafft out."
        print(result)
        if ctx.Switch_deleteIntermediateFiles == "L":
            error_resHtmlMaker(ctx.eachDirAddress, ctx.keyNode, ctx.queryID, result)
        if ctx.Switch_deleteIntermediateFiles == "D":
            deleteFiles(ctx.eachDirAddress)
        exit()

    print("\n\n##### 1st tree: TRIMAL 1st round ######\n\n")
    trimLine1 = (
        f"{ctx.toolAddress}trimal "
        f"-out {ctx.eachDirAddress}042_AA.fas.trm "
        f"-htmlout {ctx.eachDirAddress}042_AA.fas.trm.html "
        f"-in {ctx.eachDirAddress}040_mafOutAA.txt "
        f"-gappyout"
    )
    subprocess.call(trimLine1, shell=True)

    print("\n\n##### 1st tree: calcilation for ShortSequence_threshold ######\n\n")
    delete_sequences_with_alignedSiteRate(
        ctx.eachDirAddress, ctx.aligned_site_rate,
        "044_overRateAA.fas", "044_overRateDNA.fas", "044_aligned_site_rate.txt"
    )

    print("\n\n##### 1st tree: 2nd round MAFFT/trimal ######\n\n")
    res_compare_numSeqs = compare_numSeqs(
        ctx.eachDirAddress, "040_mafOutAA.txt", "044_overRateAA.fas"
    )
    if res_compare_numSeqs == "Equal":
        print("\n\n##### 2nd-round MAFFT/trimal skip and just copy files ######\n\n")
        shutil.copy(
            f"{ctx.eachDirAddress}040_mafOutAA.txt",
            f"{ctx.eachDirAddress}050_mafOutAA.txt",
        )
        shutil.copy(
            f"{ctx.eachDirAddress}042_AA.fas.trm",
            f"{ctx.eachDirAddress}052_AA.fas.trm",
        )
        shutil.copy(
            f"{ctx.eachDirAddress}042_AA.fas.trm.html",
            f"{ctx.eachDirAddress}052_AA.fas.trm.html",
        )
    else:
        print("\n\n##### 1st tree: MAFFT 2nd round ######\n\n")
        maffLine2 = (
            f"{ctx.toolAddress}mafft "
            f"{ctx.eachDirAddress}044_overRateAA.fas "
            f"> {ctx.eachDirAddress}050_mafOutAA.txt"
        )
        #subprocess.call(maffLine2, shell=True)
        run_or_make_summary_and_exit(ctx, maffLine2, "ERROR: mafft (1st round) failed")

        print("\n\n##### 1st tree: TRIMAL 2nd round ######\n\n")
        trimLine2 = (
            f"{ctx.toolAddress}trimal "
            f"-out {ctx.eachDirAddress}052_AA.fas.trm "
            f"-htmlout {ctx.eachDirAddress}052_AA.fas.trm.html "
            f"-in {ctx.eachDirAddress}050_mafOutAA.txt "
            f"-gappyout"
        )
        subprocess.call(trimLine2, shell=True)

    if ctx.dataset in ("Exclude3rd", "Include3rd"):
        print("\n\n##### 1st tree: PAL2NAL ######\n\n")
        pal2nalLine = (
            f"{ctx.toolAddress}pal2nal.pl "
            f"{ctx.eachDirAddress}050_mafOutAA.txt "
            f"{ctx.eachDirAddress}044_overRateDNA.fas "
            f"-output fasta > {ctx.eachDirAddress}054_p2nOutcDNAfas.txt"
        )
        #subprocess.call(pal2nalLine, shell=True)
        run_or_make_summary_and_exit(ctx, pal2nalLine, "ERROR: pal2nal failed")

        path_p2n = f"{ctx.eachDirAddress}054_p2nOutcDNAfas.txt"
        if os.path.getsize(path_p2n) == 0:
            makeSummary(
                ctx.SpeciesTree, ctx.taxonSamplingList, ctx.mode, ctx.dataset,
                ctx.keyNode, ctx.startTime, ctx.eachDirAddress, ctx.eachDirAddress_e1stre, ctx.BSthreshold,
                ctx.num_rootSequences, ctx.allNodes_speciesTree,
                ctx.childSpeciesNodes_AllGroup,
                aligned_site_rate=ctx.aligned_site_rate,
                outfile_summary="100_analysisSummary.txt",
            )
            result = "No pal2nal out."
            print("result", result)
            if ctx.Switch_deleteIntermediateFiles == "L":
                error_resHtmlMaker(ctx.eachDirAddress, ctx.keyNode, ctx.queryID, result)
            else:
                deleteFiles(ctx.eachDirAddress)
            exit()

    if ctx.dataset in ("Exclude3rd", "Include3rd"):
        trimaledv141_FileMakerDNA(
            ctx.eachDirAddress, "054_p2nOutcDNAfas.txt", "052_AA.fas.trm.html",
            outfile="080_trimedCDNAPhy.txt"
        )
    fas2phy(
        ctx.eachDirAddress, fastaFileName="052_AA.fas.trm",
        outPhyFileName="080_trimedAAPhy.txt"
    )

def _e_run_1st_tree_search_notung_summary(ctx):

    print("\n\n##### 1st tree: (tree search) ######")
    outgroup1 = outGroupSelect(ctx.eachDirAddress, "080_trimedAAPhy.txt")
    
    out_njbs_1st = ctx.eachDirAddress + "085_NJBS1st.txt"
    # ログは dataset によって元の命名を維持
    if ctx.dataset in ("Exclude3rd", "Include3rd"):
        log_file = "085_log.txt"
    else:
        log_file = "085_fastmelog.txt"
    log_njbs_1st = ctx.eachDirAddress + log_file
    
    if ctx.dataset == "Exclude3rd":
        print("##### APE\n\n")
        print("The 1st gene tree is estimated by excluding 3rd codon positions.")
        phyCodonToBlock(ctx.eachDirAddress, "080_trimedCDNAPhy.txt", 2, outfile="082_trimedBlockExc3rdPhy.txt")
        NJBSline1 = (
            f"{ctx.toolAddress}Rscript {ctx.scriptAddress}NJBS.R "
            f"{ctx.eachDirAddress}082_trimedBlockExc3rdPhy.txt "
            f"{outgroup1} "
            f"{out_njbs_1st} "
            f"> {log_njbs_1st} 2>&1"
        )
    elif ctx.dataset == "Include3rd":
        print("##### APE\n\n")
        print("The 1st gene tree is estimated by including 3rd codon positions.")
        phyCodonToBlock(ctx.eachDirAddress, "080_trimedCDNAPhy.txt", 3, outfile="082_trimedBlockInc3rdPhy.txt")
        NJBSline1 = (
            f"{ctx.toolAddress}Rscript {ctx.scriptAddress}NJBS.R "
            f"{ctx.eachDirAddress}082_trimedBlockInc3rdPhy.txt "
            f"{outgroup1} "
            f"{out_njbs_1st} "
            f"> {log_njbs_1st} 2>&1"
        )
    else:
        print("##### fastme\n\n")
        print("The 1st gene tree is estimated using amino acid sequences.")
        phy2fastmePhy(ctx.eachDirAddress, phyFileName="080_trimedAAPhy.txt", outFastmePhyFileName="082_trimedAAFastmePhy.txt")
        NJBSline1 = (
            f"{ctx.toolAddress}fastme -i {ctx.eachDirAddress}082_trimedAAFastmePhy.txt "
            f"--protein=WAG -m NJ -b 100 -T 2 -v 3 "
            f"-o {out_njbs_1st} "
            f"> {log_njbs_1st} 2>&1"
        )
    subprocess.call(NJBSline1, shell=True)
    if (not os.path.isfile(out_njbs_1st)) or (os.path.getsize(out_njbs_1st) == 0):
        print("[ERROR] 1st-tree search failed to create:", out_njbs_1st)
        print("[ERROR] See log:", log_njbs_1st)
        sys.exit(1)

    # 1st tree 結果の存在確認
    if not os.path.isfile(ctx.eachDirAddress + "085_NJBS1st.txt"):
        result = "No 1st tree. "
        makeSummary(
            ctx.SpeciesTree, ctx.taxonSamplingList, ctx.mode, ctx.dataset,
            ctx.keyNode, ctx.startTime, ctx.eachDirAddress, ctx.eachDirAddress_e1stre, ctx.BSthreshold,
            ctx.num_rootSequences, ctx.allNodes_speciesTree,
            ctx.childSpeciesNodes_AllGroup,
            aligned_site_rate=ctx.aligned_site_rate,
            outfile_summary="100_analysisSummary.txt",
        )
        if ctx.Switch_deleteIntermediateFiles == "L":
            error_resHtmlMaker(ctx.eachDirAddress, ctx.keyNode, ctx.queryID, result)
        if ctx.Switch_deleteIntermediateFiles == "D":
            deleteFiles(ctx.eachDirAddress)
        exit()

    nj_path = os.path.join(ctx.eachDirAddress, "085_NJBS1st.txt")
    try:
        with open(nj_path, "r") as NJtreeFile:
            NJtreeFileCont = list(NJtreeFile)
    except OSError:
        NJtreeFileCont = []
    if not NJtreeFileCont:
        result = "No 1st tree. "
        makeSummary(
            ctx.SpeciesTree, ctx.taxonSamplingList, ctx.mode, ctx.dataset,
            ctx.keyNode, ctx.startTime, ctx.eachDirAddress, ctx.eachDirAddress_e1stre, ctx.BSthreshold,
            ctx.num_rootSequences, ctx.allNodes_speciesTree,
            ctx.childSpeciesNodes_AllGroup,
            aligned_site_rate=ctx.aligned_site_rate,
            outfile_summary="100_analysisSummary.txt",
        )
        if ctx.Switch_deleteIntermediateFiles == "L":
            error_resHtmlMaker(ctx.eachDirAddress, ctx.keyNode, ctx.queryID, result)
        if ctx.Switch_deleteIntermediateFiles == "D":
            deleteFiles(ctx.eachDirAddress)
        sys.exit(1)

    print("\n\n##### 1st tree: NOTUNG ######\n\n")
    NOTUNG1stLine = (
        f"java -jar {ctx.toolAddress}Notung.jar "
        f"-s {ctx.eachDirAddress}000_speciesTree_topLeft.txt "
        f"-g {ctx.eachDirAddress}085_NJBS1st.txt --outputdir {ctx.eachDirAddress} "
        f"--rearrange --threshold {ctx.BSthreshold} --speciestag prefix "
        f"--maxtrees 5 --nolosses --treeoutput nhx > {ctx.eachDirAddress}085_NOTUNGlog.txt"
    )
    print("NOTUNG1stLine:", NOTUNG1stLine)
    subprocess.call(NOTUNG1stLine, shell=True)

    if not os.path.isfile(ctx.eachDirAddress + "085_NJBS1st.txt.rearrange.0"):
        print("Error in NOTUNG: Cannot compare the NJ and species tree.")
        print("Check your species tree.")
        exit()

    print("\n\n##### 1st tree: Making summary ######\n\n")
    makeSummary(
        ctx.SpeciesTree, ctx.taxonSamplingList, ctx.mode, ctx.dataset,
        ctx.keyNode, ctx.startTime, ctx.eachDirAddress, ctx.eachDirAddress_e1stre, ctx.BSthreshold, ctx.BSthreshold_4_ReE1st,
        ctx.num_rootSequences, ctx.allNodes_speciesTree,
        ctx.childSpeciesNodes_AllGroup,
        aligned_site_rate=ctx.aligned_site_rate,
        outfile_summary="100_analysisSummary.txt",
    )

    if ctx.Switch_deleteIntermediateFiles == "L":
        print("\n\n##### 1st tree: APE (tree draw) -2 ######\n\n")
        if ctx.mode == "E":
            thick_branch_species_name_line = "Orthogroup"
            root_species_name_line = "Rooting_4_2ndTree"
        elif ctx.mode == "E1st":
            thick_branch_species_name_line = "Rooting_4_1stTree"
            root_species_name_line = "Rooting_4_1stTree"
        else:
            print("Error. Check mode name: ", ctx.mode)
            exit()
        treePlotR_1st2 = (
            f"{ctx.toolAddress}Rscript {ctx.scriptAddress}treePlot.R "
            f"{ctx.eachDirAddress}100_analysisSummary.txt "
            f"1st_gene_tree_newick "
            f"1st_rearranged_gene_tree_newick "
            f"{thick_branch_species_name_line} "
            f"{root_species_name_line} "
            f"{ctx.eachDirAddress}115_1st > {ctx.eachDirAddress}115_logTreePlotB.txt"
        )
        print("treePlotR_1st2: ", treePlotR_1st2)
        subprocess.call(treePlotR_1st2, shell=True)

        if ctx.mode == "E1st":
            make_resHtml2(
                ctx.queryID,
                ctx.mode,
                ctx.BSthreshold_4_ReE1st,
                ctx.eachDirAddress,
                resHTMLlines_2steps
            )

def _e_run_2nd_tree_pipeline(ctx):

    # 4) 2nd tree
    print("\n\n##### 2nd tree ######\n\n")

    resDict_1st = readRes_dict(ctx.eachDirAddress + "100_analysisSummary.txt")
    ctx.resDict_1st = resDict_1st

    if ctx.resDict_1st[">BS_of_orthogroupBasalNode"][0].startswith("noOrthogroup_"):
        print("### 1st tree: noOrthogroup_, then stop here")
        result = ctx.resDict_1st[">BS_of_orthogroupBasalNode"][0]
        if ctx.Switch_deleteIntermediateFiles == "L":
            error_resHtmlMaker(ctx.eachDirAddress, ctx.keyNode, ctx.queryID, result)
        if ctx.Switch_deleteIntermediateFiles == "D":
            deleteFiles(ctx.eachDirAddress)
        exit()

    if len(ctx.resDict_1st[">Orthogroup"]) < 4:
        result = "Less than 4 orthogroup members."
        if ctx.Switch_deleteIntermediateFiles == "L":
            error_resHtmlMaker(ctx.eachDirAddress, ctx.keyNode, ctx.queryID, result)
        if ctx.Switch_deleteIntermediateFiles == "D":
            deleteFiles(ctx.eachDirAddress)
        exit()

    rootLeaves = selectRootSp4secondTreeSearch(ctx.eachDirAddress, ctx.keyNode, ctx.num_rootSequences)

    if ctx.dataset in ("Exclude3rd", "Include3rd"):
        make_2ndanalysis_seqFile(ctx.eachDirAddress, rootLeaves, "054_p2nOutcDNAfas.txt", ctx.resDict_1st, outfile="150_noGapCDNAfas.txt")
        cDNAfas2noGapAAFasFile(ctx.eachDirAddress, "150_noGapCDNAfas.txt", outfile="150_noGapAA.txt")
    else:
        make_2ndanalysis_seqFile(ctx.eachDirAddress, rootLeaves, "050_mafOutAA.txt", ctx.resDict_1st, outfile="150_noGapAA.txt")

    print("\n\n##### 2nd tree: MAFFT ######\n\n")
    maffLine2 = (
        f"mafft {ctx.eachDirAddress}150_noGapAA.txt "
        f"> {ctx.eachDirAddress}160_mafOut.txt"
    )
    #subprocess.call(maffLine2, shell=True)
    run_or_make_summary_and_exit(ctx, maffLine2, "ERROR: mafft (1st round) failed")

    print("\n\n##### 2nd tree: TRIMAL ######\n\n")
    trimLine2 = (
        f"{ctx.toolAddress}trimal "
        f"-out {ctx.eachDirAddress}170_trimedAAOutFas.txt "
        f"-htmlout {ctx.eachDirAddress}170_aln_prot.html "
        f"-in {ctx.eachDirAddress}160_mafOut.txt "
        f"-gappyout"
    )
    subprocess.call(trimLine2, shell=True)

    fas2phy(ctx.eachDirAddress, fastaFileName="160_mafOut.txt", outPhyFileName="190_aln_prot.txt")
    fas2phy(ctx.eachDirAddress, "170_trimedAAOutFas.txt", "200_trimedAAPhy.txt")

    if ctx.dataset in ("Exclude3rd", "Include3rd"):
        print("\n\n##### 2nd tree: PAL2NAL ######\n\n")
        pal2nalLine = (
            f"{ctx.toolAddress}pal2nal.pl "
            f"{ctx.eachDirAddress}160_mafOut.txt "
            f"{ctx.eachDirAddress}150_noGapCDNAfas.txt "
            f"-output fasta > {ctx.eachDirAddress}180_aln_nucl_fas.txt"
        )
        #subprocess.call(pal2nalLine, shell=True)
        run_or_make_summary_and_exit(ctx, pal2nalLine, "ERROR: pal2nal failed")
        


        fas2phy(ctx.eachDirAddress, fastaFileName="180_aln_nucl_fas.txt", outPhyFileName="190_aln_nucl.txt")

        trimaledv141_FileMakerDNA(ctx.eachDirAddress, "180_aln_nucl_fas.txt", "170_aln_prot.html", outfile="200_trimedCDNAPhy.txt")

        phyCodonToBlock(ctx.eachDirAddress, "200_trimedCDNAPhy.txt", 2, outfile="210_trimedBlockExc3rdPhy.txt")
        phyCodonToBlock(ctx.eachDirAddress, "200_trimedCDNAPhy.txt", 3, outfile="210_trimedBlockInc3rdPhy.txt")

    # 2nd tree 検索
    print("\n\n##### 2nd tree: APE (tree search) ######\n\n")
    resDict_OG = readRes_dict(ctx.eachDirAddress + "100_analysisSummary.txt")
    outGroup = resDict_OG[">Rooting_4_2ndTree"][0]
    print("treeSearchMethod", ctx.treeSearchMethod)
    if ctx.treeSearchMethod == "ML":
        make_raxmlPartitionFile(ctx.eachDirAddress, outPartFile="220_raxmlPart.txt")
        raxmlLine = (
            f"{ctx.scriptAddress}raxmlHPC-PTHREADS-SSE3 "
            f"-f a -x 12345 -p 12345 -# 100 "
            f"-m GTRCAT "
            f"-s {ctx.eachDirAddress}210_trimedBlockExc3rdPhy.txt "
            f"-w {ctx.eachDirAddress} "
            f"-q {ctx.eachDirAddress}220_raxmlPart.txt "
            f"-o {outGroup} "
            f"-n txt "
            f"-T 2"
        )
        subprocess.call(raxmlLine, shell=True)
        moveRAxMLfiles(ctx.eachDirAddress, outfile="230_2ndtree.txt")
    else:
        if ctx.dataset == "Exclude3rd":
            print("The 2nd gene tree is estimated by excluding 3rd codon positions.")
            NJBSline2 = (
                f"{ctx.toolAddress}Rscript {ctx.scriptAddress}NJBS.R "
                f"{ctx.eachDirAddress}210_trimedBlockExc3rdPhy.txt "
                f"{outGroup} "
                f"{ctx.eachDirAddress}230_2ndtree.txt "
                f"> {ctx.eachDirAddress}230_log.txt"
            )
        elif ctx.dataset == "Include3rd":
            print("The 2nd gene tree is estimated by including 3rd codon positions.")
            NJBSline2 = (
                f"{ctx.toolAddress}Rscript {ctx.scriptAddress}NJBS.R "
                f"{ctx.eachDirAddress}210_trimedBlockInc3rdPhy.txt "
                f"{outGroup} "
                f"{ctx.eachDirAddress}230_2ndtree.txt "
                f"> {ctx.eachDirAddress}230_log.txt"
            )
        else:
            print("The 2nd gene tree is estimated using amino acid sequences.")
            phy2fastmePhy(ctx.eachDirAddress, phyFileName="200_trimedAAPhy.txt", outFastmePhyFileName="202_trimedAAFastmePhy.txt")
            NJBSline2 = (
                f"{ctx.toolAddress}fastme "
                f"-i {ctx.eachDirAddress}202_trimedAAFastmePhy.txt "
                f"--protein=WAG -m NJ -b 100 -v 3 "
                f"-o {ctx.eachDirAddress}230_2ndtree.txt "
                f"> {ctx.eachDirAddress}230_log.txt"
            )
        subprocess.call(NJBSline2, shell=True)

    if not os.path.isfile(ctx.eachDirAddress + "230_2ndtree.txt"):
        print("Error: Cannot estimate 2nd tree.")
        if ctx.Switch_deleteIntermediateFiles == "D":
            deleteFiles(ctx.eachDirAddress)
        exit()

    Tree2ndFilePath = ctx.eachDirAddress + "230_2ndtree.txt"
    if os.path.getsize(Tree2ndFilePath) == 0:
        result = "2nd Tree Error: Cannot estimate the 2nd NJ tree."
        print("result", result)
        make_resHtml2(
            ctx.queryID,
            ctx.mode,
            ctx.BSthreshold_4_ReE1st,
            ctx.eachDirAddress,
            resHTMLlines_2steps
        )
        if ctx.Switch_deleteIntermediateFiles == "D":
            deleteFiles(ctx.eachDirAddress)
        exit()

    rootBS100R = (
        f"{ctx.toolAddress}Rscript {ctx.scriptAddress}rootBS100.R "
        f"{ctx.eachDirAddress}230_2ndtree.txt "
        f"{ctx.eachDirAddress}230_2ndtreeRootBS100.txt"
    )
    subprocess.call(rootBS100R, shell=True)

    print("\n\n##### 2nd tree: NOTUNG ######\n\n")
    NOTUNG2ndLine = (
        f"java -jar {ctx.toolAddress}Notung.jar "
        f"-s {ctx.eachDirAddress}000_speciesTree_topLeft.txt "
        f"-g {ctx.eachDirAddress}230_2ndtreeRootBS100.txt "
        f"--outputdir {ctx.eachDirAddress} "
        f"--rearrange "
        f"--threshold {ctx.BSthreshold} "
        f"--speciestag prefix "
        f"--maxtrees 5 "
        f"--nolosses "
        f"--treeoutput nhx "
        f"> {ctx.eachDirAddress}230_NOTUNGlog.txt"
    )
    subprocess.call(NOTUNG2ndLine, shell=True)


    add_makeSummary(
        ctx.eachDirAddress, ctx.treeSearchMethod, ctx.keyNode,
        ctx.allNodes_speciesTree, ctx.childSpeciesNodes_focalGroup,
        ctx.startTime, outfile_summary2="100_analysisSummary.txt",
    )

def _e_finalize(ctx):
    if ctx.Switch_deleteIntermediateFiles == "L":
        print("\n\n##### 2nd tree: APE (tree draw) ######\n\n")
        treePlotR_2nd = (
            f"{ctx.toolAddress}Rscript {ctx.scriptAddress}treePlot.R "
            f"{ctx.eachDirAddress}100_analysisSummary.txt "
            f"2nd_gene_tree_newick "
            f"2nd_rearranged_gene_tree_newick "
            f"Rooting_4_2ndTree "
            f"Rooting_4_2ndTree "
            f"{ctx.eachDirAddress}240_2nd > {ctx.eachDirAddress}240_logTreePlotB.txt"
        )
        print("treePlotR_2nd: ", treePlotR_2nd)
        subprocess.call(treePlotR_2nd, shell=True)
        
        make_resHtml2(
            ctx.queryID,
            ctx.mode,
            ctx.BSthreshold_4_ReE1st,
            ctx.eachDirAddress,
            resHTMLlines_2steps
        )
    phy2fas(ctx.eachDirAddress, infile_phy="190_aln_prot.txt", outfile_fas="190_aln_prot_fas.txt")
    copy_alignment_orthogroup(ctx.eachDirAddress, ctx.dataset, ctx.queryID, ctx.alignment_orthogroups)

    if ctx.Switch_deleteIntermediateFiles == "D":
        deleteFiles(ctx.eachDirAddress)

    elapsed = time.time() - ctx.startTime
    print("\n\n#####")
    print(f"Mode {ctx.mode} analysis has completed successfully. Elapsed time: {elapsed:.1f} seconds.\n\n")
    exit()


### File checking
def check_scriptsDirectory(scriptAddress):
    """
    scripts ディレクトリに必要な R スクリプトがあるかを確認する。
    - 必須（完全一致）: NJBS.R, ladderizeTree.R, treePlot.R, rootBS100.R
    見つからない場合は、明確なエラーメッセージを出して終了。
    """
    import os, sys
    if not os.path.exists(scriptAddress):
        print("Error. Your >scripts is not found:")
        print(scriptAddress)
        sys.exit(1)

    required_exact = ["NJBS.R", "ladderizeTree.R", "treePlot.R", "rootBS100.R"]

    missing = []
    for fname in required_exact:
        if not os.path.isfile(os.path.join(scriptAddress, fname)):
            missing.append(fname)

    if missing:
        print("Error: Required R scripts are missing in your scripts directory:", scriptAddress)
        for m in missing:
            print("-", m)
        print("Please place the files above under the scripts directory specified by >scripts in control.txt.")
        sys.exit(1)

def check_toolsDirectory(mode, toolAddress, dataset):
    dependencies = os.listdir(path=toolAddress)
    #print("dependencies", dependencies)
    #exit()
    if mode == "E":
        if not "blastp" in dependencies:
            print("Error: Cannot find blastp in your tools directory:", toolAddress)
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
        if not "makeblastdb" in dependencies:
            print("Error: Cannot find makeblastdb in your tools directory:", toolAddress)
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
        if not "mafft" in dependencies:
            print("Error: Cannot find mafft in your tools directory:", toolAddress)
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
        if "trimal" not in dependencies:
            print("Error: Cannot find trimal in your tools directory:", toolAddress)
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
        if dataset == "AminoAcid":
            if not "fastme" in dependencies:
                print("Error: Cannot find fastme in your tools directory:", toolAddress)
                print(f"Dataset {dataset} was selected in the control.txt file.")
                print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
                sys.exit()
        if not "pal2nal.pl" in dependencies:
            print("Error: Cannot find pal2nal.pl in your tools directory:", toolAddress)
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
        if not "Notung.jar" in dependencies:
            print("Error: Cannot find Notung.jar in your tools directory:", toolAddress)
            print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
            exit()
    if not "Rscript" in dependencies:
        print("Error: Cannot find Rscript in your tools directory:", toolAddress)
        print("See the tutorial from https://github.com/jun-inoue/ORTHOSCOPE_STAR")
        exit()


def check_mode(mode, queryID):
    #print("### check_mode() ###")
    #print(mode)
    #exit()
    # Mode E/E1st/ReE1st/D は「引数は gene ID（.txt ではない）」を想定
    if mode in ("E", "E1st", "ReE1st", "D", "D1st", "ReD1st") and queryID.endswith(".txt"):
        print("Error in your control.txt file.")
        print(f"You selected mode {mode}. Mode {mode} needs gene ID as the argument.")
        print("Your argument was:", queryID)
        sys.exit(1)

    # Mode S/S1st は「引数は .txt ファイル」を想定
    if mode in ("S", "S1st") and not queryID.endswith(".txt"):
        print("Error in your control.txt file.")
        print(f"You selected mode {mode}. Mode {mode} needs .txt file as the argument.")
        print("Your argument was:", queryID)
        sys.exit(1)

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
    return re.sub(r" ", "", resDict_SR[keyword][0])


def check_pickup_taxonSampling(dbAddress, lines_taxonSampling):
    dbLines = []
    taxonSamplingList = []
    for line in lines_taxonSampling:
        lineTMP = re.sub(r" +", " ", line)
        lineTMP = re.sub(r"^ +", "", lineTMP)
        lineTMP = re.sub(r" +$", "", lineTMP)

        line_separated = lineTMP.split(" ")
        if len (line_separated) != 3:
            print("Error in your TaxonSampling in the control file.")
            print("Each line should consists of three part separated by spaces.")
            print("Check the following line:")
            print(re.sub(r" +", " ", line))
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


def check_presense_of_databases(dbAddress, dbLines):
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
    name_lines = [
        ">QuerySpecies",
        ">Mode",
        ">TaxonSampling",
        ">SpeciesTree",
        ">Num_rootSequences",
        ">KeyNode",
        ">BLAST_Evalue_threshold_for_reported_sequences",
        ">Number_of_hits_to_report_per_genome",
        ">ShortSequence_threshold",
        ">Dataset",
        ">BSthreshold",
        ">BSthreshold_4_ReE1st",
        ">Outdir",
        ">Outdir_4_ReE1st",
        ">Database",
        ">tools",
        ">scripts"
    ]
    for name_line in name_lines:
        if not name_line in resDict_SR.keys():
            print(f"Error in your control file. Nameline {name_line} is not found.")
            print("From v1.2.4, >tools and >scripts needs to be defined in the control.txt file.")
            exit()


def read_controlFile():
    print("### read_controlFile() ###")
    resDict_SR = readRes_dict("control.txt")        
    check_controlFile(resDict_SR)

    ## the others
    SpeciesTree = check_pickup_parameter(resDict_SR, "SpeciesTree")
    blastEvalue = check_pickup_parameter(resDict_SR, "BLAST_Evalue_threshold_for_reported_sequences")
    if re.search(r"[^\de-]", blastEvalue):
        print("Error. >BLAST_Evalue_threshold_for_reported_sequences should be floating point expression (e.g., 1e-3) or just 1.")
        print("Your >BLAST_Evalue_threshold_for_reported_sequences is", blastEvalue)
        exit()

    Number_of_hits_to_report_per_genome = check_pickup_parameter(resDict_SR, "Number_of_hits_to_report_per_genome")
    if re.search(r"[^\d]", Number_of_hits_to_report_per_genome):
        print("Error. >Number_of_hits_to_report_per_genome should be integer.")
        print("Your >Num_rootSequences is", Number_of_hits_to_report_per_genome)
        exit()

    aligned_site_rate = check_pickup_parameter(resDict_SR, "ShortSequence_threshold")
    if re.search(r"[^\d\.]", aligned_site_rate):
        print("Error. >ShortSequence_threshold should be floating less than 1.")
        print("Your >Num_rootSequences is", aligned_site_rate)
        exit()

    dataset = check_pickup_parameter(resDict_SR, "Dataset")
    if dataset == "Exclude3rd" or dataset == "Include3rd" or dataset == "AminoAcid":
        pass
    else:
        print("Error. >Dataset should be Exclude3rd, Include3rd, or AminoAcid.")
        print("Your >Dataset is", dataset)
        exit()

    # ---- BSthreshold ----------------------------------------------------
    BSthreshold = check_pickup_parameter(resDict_SR, "BSthreshold")
    BSthreshold = BSthreshold.strip()
    if BSthreshold.isdigit():
        val = int(BSthreshold)
        if not (0 <= val <= 100):
            print("Error. >BSthreshold should be an integer between 0 and 100.")
            print("Your >BSthreshold is", BSthreshold)
            exit()
    else:
        if BSthreshold not in ("reconcile", "transfer"):
            print("Error. >BSthreshold should be 0–100 or 'reconcile' or 'transfer'.")
            print("Your >BSthreshold is", BSthreshold)
            exit()
    
    # ---- BSthreshold_4_ReE1st ------------------------------------------
    BSthreshold_4_ReE1st = check_pickup_parameter(resDict_SR, "BSthreshold_4_ReE1st")
    BSthreshold_4_ReE1st = BSthreshold_4_ReE1st.strip()
    if BSthreshold_4_ReE1st.isdigit():
        val = int(BSthreshold_4_ReE1st)
        if not (0 <= val <= 100):
            print("Error. >BSthreshold_4_ReE1st should be an integer between 0 and 100.")
            print("Your >BSthreshold_4_ReE1st is", BSthreshold_4_ReE1st)
            exit()
    else:
        if BSthreshold_4_ReE1st not in ("reconcile", "transfer"):
            print("Error. >BSthreshold_4_ReE1st should be 0–100 or 'reconcile' or 'transfer'.")
            print("Your >BSthreshold_4_ReE1st is", BSthreshold_4_ReE1st)
            exit()

    #treeSearchMethod = check_pickup_parameter(resDict_SR, "TreeSearchMethod")
    treeSearchMethod = "NJ"
    num_rootSequences = check_pickup_parameter(resDict_SR, "Num_rootSequences")
    if re.search(r"[^\d]", num_rootSequences):
        print("Error. >Num_rootSequences should be integer. 1 is favorable to estimate sister groups.")
        print("Your >Num_rootSequences is", num_rootSequences)
        exit()

    keyNode = check_pickup_parameter(resDict_SR, "KeyNode")
    #name_querySpeciesNode = check_pickup_parameter(resDict_SR, "QuerySpeciesGroup")
    name_querySpecies = check_pickup_parameter(resDict_SR, "QuerySpecies")
    speciesWithGeneFunction = check_pickup_parameter(resDict_SR, "SpeciesWithGeneFunction")
    outdir = check_pickup_parameter(resDict_SR, "Outdir")
    outdir_4_ReE1st = check_pickup_parameter(resDict_SR, "Outdir_4_ReE1st")
    alignment_orthogroups = check_pickup_parameter(resDict_SR, "Alignment_orthogroups")
    mode = check_pickup_parameter(resDict_SR, "Mode")

    allowed_modes = ("E", "E1st", "ReE1st", "D", "D1st", "ReD1st", "S", "S1st", "ReS1st")
    if mode in allowed_modes:
        pass
    else:
        print("Error. Check your >Mode.")
        print(f">Mode should be one of: {', '.join(allowed_modes)}.")
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
    #print("dbAddress", dbAddress)
    if not os.path.exists(dbAddress):
        print("Error. Your >Database is not found:")
        print(dbAddress)
        exit()
    dbAddress = dbAddress + "/"

    toolAddress = check_pickup_parameter(resDict_SR, "tools")
    #print("dbAddress", dbAddress)
    if not os.path.exists(dbAddress):
        print("Error. Your >tools is not found:")
        print(dbAddress)
        exit()
    toolAddress = toolAddress + "/"

    scriptAddress = check_pickup_parameter(resDict_SR, "scripts")
    #print("dbAddress", dbAddress)
    if not os.path.exists(dbAddress):
        print("Error. Your >scripts is not found:")
        print(dbAddress)
        exit()
    scriptAddress = scriptAddress + "/"

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

    return (
        dbLinesTMP,
        taxonSamplingListTMP,
        SpeciesTree,
        blastEvalue,
        Number_of_hits_to_report_per_genome,
        aligned_site_rate,
        dataset,
        BSthreshold,
        BSthreshold_4_ReE1st,
        treeSearchMethod,
        num_rootSequences,
        keyNode,
        name_querySpecies,
        queryDatabase,
        dbAddress,
        toolAddress,
        scriptAddress,
        outdir,
        outdir_4_ReE1st,
        alignment_orthogroups,
        mode,
        Switch_deleteIntermediateFiles,
        speciesWithGeneFunction,
    )

def reorder_dbLines(dbLinesTMP, taxonSamplingListTMP, name_querySpecies, dirAddress_FN, treeFileName):
    #print("#### reorder_dbLines ####")
    #print("name_querySpecies", name_querySpecies)
    leaves = collect_leaves_InOrderFrom_bothNHXnewick(dirAddress_FN, treeFileName)

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


def dirFileMake(mode, outdir, eachDirAddress, alignment_orthogroups):
    #print("### dirFileMake ###")
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if not os.path.exists(eachDirAddress):
        os.mkdir(eachDirAddress)

    if mode == "E":
        if not os.path.exists(alignment_orthogroups):
            os.mkdir(alignment_orthogroups)


def make_querySeqFile(queryDatabase, name_querySpecies, dbAddress, eachDirAddress, queryID):
    #print("### make_querySeqFile ###")
    #print("queryDatabase", queryDatabase)
    #print("name_querySpecies", name_querySpecies)
    #print("queryID", queryID)
    #exit()

    recs_queryDB = readFasta_dict(dbAddress, queryDatabase)

    nameline_query = ""
    sequence_query = ""
    for nameline_DB, sequence_DB in recs_queryDB.items():
        #print("nameline_DB", nameline_DB)
        if re.search(r">" + queryID + " ", nameline_DB) or re.search(r">" + queryID + "_", nameline_DB) or re.search(r">" + queryID + "$", nameline_DB):
            nameline_query = nameline_DB
            sequence_query = sequence_DB
            break
    if not nameline_query:
        print("Error in your control file:")
        print(queryID, "is not found in", dbAddress + queryDatabase)
        exit()
    
    fs = open(eachDirAddress + "000_aaSeq_assigned_by_ID.txt", "w")
    nameline_query = re.sub(r" .*$", "", nameline_query)
    nameline_query = re.sub(r">", ">" + name_querySpecies + "_", nameline_query)
    fs.write(nameline_query + "\n")
    fs.write(sequence_query + "\n")
    fs.close()


def makeblastdb_database(dbAddress, dbLines, toolAddress):
    # names_dbFile = files = os.listdir(dbAddress)
    extensions = ["phr", "pin", "psq"]
    for dbLine in dbLines:
        flag = 0
        for extension in extensions:
            if os.path.isfile(dbAddress + dbLine[1] + "." + extension):
                flag += 1
        if flag < 3:  # ← HTMLエスケープではなく < を使う
            print("#### Running makeblastdb.")
            # 末尾にスラッシュが入っている toolAddress を前提に連結
            line_makeblastdb = (
                f"{toolAddress}"
                f"makeblastdb -dbtype prot -in "
                f"{dbAddress}{dbLine[1]}"
            )
            subprocess.call(line_makeblastdb, shell=True)

### File checking End
##############################################


##############################################
### cDNA and AA file make Start
def check_uploaded_file_as_fasta_format(eachDirAddress):
    #print("### check_uploaded_file_as_fasta_format() ###")
    f = open(eachDirAddress + "000_aaSeq_assigned_by_ID.txt")
    lines = list(f)
    f.close()
    if not lines:
        print('Error: Cannot find your sequence file in the specified directory.\n')
        exit()
    if not lines[0].startswith(">"):
        print('Error: Check your sequence file. In fasta format, name line starts with ">"\n')
        exit()
    #recs_uploaded = readFasta_dict(eachDirAddress , "000_aaSeq_assigned_by_ID.txt")
    #for name, seq in recs_uploaded.items():
    #    if re.search(r"[^ATGCNX ]", seq):
    #        print ('Error: Check your sequence. Sequences should be consist of A,T,G,C,N,X.')
    #        print ("in ", name)
    #        exit()


def ckeck_cDNAsequence(recsFN):
    for name, seq in recsFN.items():
        if re.search(r"[^ATGCNXatgcnx ]", seq):
            print ("Error: Please check sequence in: ", name)
            print ("The sequence should be A,T,G,C,N,X. Analysis stopped.")
            exit()


def readPhy_dict(eachDirAddress, phyFileName):
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
        if re.search(r"^#", Line) or re.search(r"^$", Line):
            continue
        Line = Line.rstrip("\n")
        if Line[0] == ">":
            if flag == 0:
                Name = Line
                flag = 1
            else:
                Name = re.sub(r" +$", "", Name)
                seqDictFN[Name] = stock
                Name = Line
                stock = []
        else:
            stock.append(Line)
    Name = re.sub(r" +$", "", Name)
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
            Name = Line
            Name = re.sub(r" +$", "", Name)
            seqDictFN[Name] = ""
        else:
            Line = Line.replace("\n", "")
            Line = Line.replace("\r", "")
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


def aaSeqMaker(eachDirAddress):
    recfn = readFasta_dict(eachDirAddress , "000_aaSeq_assigned_by_ID.txt")

    fa = open(eachDirAddress + "000_translated_cds_assigned_by_ID.txt","w")

    # >Human_ENSP00000259365 gene:ENSG00000136842 transcript:ENST00000259365 gene_biotype:protein_coding
    for name, seq in recfn.items():
        newName = re.sub(r"\|", "", name)
        fa.write(newName + "\n")
        #fa.write(translation(seq) + "\n")
        fa.write(seq + "\n")
    fa.close()


### cDNA and AA file make End
##############################################


###############################################################
##################### Tree manipulation START
def get_species_tree_info(SpeciesTree, keyNode, name_querySpecies):
    """
    種系統樹（SpeciesTree）に基づくノード情報を収集して返す純粋関数。
    ctx は参照・更新しない。呼び出し側で戻り値を ctx にセットする。

    Parameters
    ----------
    SpeciesTree : str
        種系統樹（本スクリプトが想定する Newick/NHX 由来表現）
    keyNode : str
        解析のキーとなるノード名
    name_querySpecies : str
        クエリ種名

    Returns
    -------
    tuple
        (
            allNodes_speciesTree,
            focalNode_speciesTree,
            speciesNodes_including_querySpecies,
            childSpeciesNodes_AllGroup,
            childSpeciesNodes_focalGroup,
        )
    """

    if not SpeciesTree or not isinstance(SpeciesTree, str):
        raise ValueError("SpeciesTree is empty or not a string.")
    if not keyNode:
        raise ValueError("keyNode is empty.")
    if not name_querySpecies:
        raise ValueError("name_querySpecies is empty.")


    # 1) 全ノードの収集
    allNodes_speciesTree = collect_nodes_from_speciesTree(SpeciesTree)

    # 2) keyNode に対応するフォーカル種ノード
    focalNode_speciesTree = identifiy_focalNode_speciesTree(
        allNodes_speciesTree, keyNode
    )

    # 3) クエリ種ノードとその祖先系列
    recs_querySpeciesNode = identify_speciesNode(
        allNodes_speciesTree, name_querySpecies
    )
    speciesNodes_including_querySpecies = collect_ancestralNodes(
        allNodes_speciesTree, recs_querySpeciesNode
    )

    # 4) keyNode 配下でクエリ種を含む子ノード（全体）
    childSpeciesNodes_AllGroup = collect_AllchildSpeciesNodes_with_querySpecies(
        allNodes_speciesTree, focalNode_speciesTree, name_querySpecies
    )

    # 5) フォーカル群の中でクエリ種を含む子ノード
    childSpeciesNodes_focalGroup = collect_childSpeciesNodes_with_querySpecies(
        allNodes_speciesTree, focalNode_speciesTree, name_querySpecies
    )

    return (
        allNodes_speciesTree,
        focalNode_speciesTree,
        speciesNodes_including_querySpecies,
        childSpeciesNodes_AllGroup,
        childSpeciesNodes_focalGroup,
    )

def collect_leaves_InOrderFrom_bothNHXnewick(eachDirAddress, treeFileName):
    #print("### collect_leaves_InOrderFrom_bothNHXnewick() ###")
    treeTMP = open(eachDirAddress + treeFileName, "r")
    tree = list(treeTMP)[0]
    treeTMP.close()
    #print("tree", tree)
    leaves = tree.split(",")
    #leaves = [re.sub(r"^\(*([^:\(\)]+).*$", r"\1", leaf) for leaf in leaves]
    for i, leaf in enumerate(leaves):
        leaves[i] = re.sub(r"^\(*([^:\(\)]+).*$", r"\1", leaf)
    #leaves = [re.sub(r"[\t\n]", "", leaf) for leaf in leaves]
    for i, leaf in enumerate(leaves):
        leaves[i] = re.sub(r"[\t\n]", "", leaf)
    #for leaf in leaves:
    #    print("leaf", leaf)
    #exit()
    return leaves


def test_querySpecies_is_in_speciesTree(taxonSamplingList, SpeciesTree):
    speciesNames_taxonSamplingList = [] 
    for speciesName_taxonSamplingList in taxonSamplingList:
        speciesName_taxonSamplingList = re.sub(r"_[^_]+$", "", speciesName_taxonSamplingList)
        speciesNames_taxonSamplingList.append(speciesName_taxonSamplingList)

    if re.search(r"\)\)", SpeciesTree) or re.search(r"\),", SpeciesTree) :
        print ("Error: In the species tree, all nodes should have node name.")
        print ("Please check you newick format using FigTree or TreeGraph_2.")
        exit()

    num_rightOpen = SpeciesTree.count("(")
    num_leftOpen  = SpeciesTree.count(")")
    num_comma     = SpeciesTree.count(",")
    if num_rightOpen != num_leftOpen:
        print ("Error: In the species tree, numbers of right and left parehtheses should be same.")
        exit()
    if num_rightOpen != num_comma or num_leftOpen != num_comma :
        print ("Error: In the species tree, all nodes should be bufuricated.")
        exit()

    speciesNames_in_speciesTree = ""
    for node in allNodes_speciesTree:
        #print("node", node)
        if re.search(r"S=" + keyNode + ":", node[2]):
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
        spName_dbLine = re.sub(r"_$","",spName_dbLine)
        if not spName_dbLine in allSpeciesNames_speciesTree:
            print ("Error: ", spName_dbLine, ", is not found in your uploaded-species tree.")
            exit()


def change_nhx_to_newick_with_NHXnodeName(tree_NHX):
  
    cladeReg = r"\)([^\[]+\[&&NHX[^\]]+\])"
    count = 0
    while re.search(cladeReg, tree_NHX):
        #print("count   :", count)
        matchA = re.search(cladeReg, tree_NHX)
        expTemp  = matchA.group(1);
        #print("expTemp:", expTemp)
        #if count == 10:
        #    exit()

        matchB = re.search(r"(\[.*\])", expTemp)
        exp  = matchB.group(1)
 
        if re.search(r"B=([\d]+)", expTemp):
            matchC = re.search(r"B=([\d]+)", expTemp)
            bs = matchC.group(1)
        else:
            bs = "r"

        exp = re.sub(r"&&NHX:", "", exp)
        exp = re.sub(r":",      "_", exp)
        exp = re.sub(r"_B=.*$", "", exp)
        
        #print("expTemp :", expTemp)
        #print("bs      :", bs)
        #print("exp     :", exp)
        #print("bs + exp:", bs + exp)

        tree_NHX = re.sub(cladeReg, ")" + bs + "_" + exp, tree_NHX, count=1)
        count += 1
        #print("tree2: ", tree_NHX)
        #print()
    
    tree_NHX = re.sub(r"\:[\d|\.E-]+", "", tree_NHX)   #Delete the branch lengths

    ## left node name
    tree_NHX = re.sub(r"\[&&NHX:[^]]+\]", "",  tree_NHX)
    tree_NHX = re.sub(r"[\[\]]",          "", tree_NHX)
    #tree_NHX = re.sub(r"\]",              "",  tree_NHX)

    return tree_NHX;


def collect_nodes_from_speciesTree(SpeciesTree):
    #print("### collect_nodes_from_speciesTree() ###")
    #print("SpeciesTree", SpeciesTree)
    #exit()
    tree_newick = SpeciesTree
    nodes = []     # 2D array for nodes
    cladeReg = r"\(([^\(\)]+)\)(\w+)"
    while re.search(cladeReg, tree_newick):
        tree_newick = tree_newick.rstrip("\n")
        match = re.search(cladeReg, tree_newick)             # Pick up the smallest and leftmost clade for the following analysis
        leavesString = match.group(1)
        exp = match.group(2)
        #leaves = set(leavesString.split(","))
        leaves = leavesString.split(",")
        tree_newick = re.sub(cladeReg, r"\1", tree_newick, count=1)  # Delete the outer parentheses from the analyzed clade
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


#def collect_nodes_from_newick(tree_newick):
#
#    nodes = []     # 2D array for nodes
#    cladeReg = "\(([^\(\)]+)\)(\w+)"
#    while re.search(cladeReg, tree_newick):
#        tree_newick = tree_newick.rstrip("\n")
#        match = re.search(cladeReg, tree_newick)             # Pick up the smallest and leftmost clade for the following analysis
#        leavesString = match.group(1)
#        exp = match.group(2)
#        #leaves = set(leavesString.split(","))
#        leaves = leavesString.split(",")
#        #print("leaves", leaves)
#        #exit()
#        tree_newick = re.sub(cladeReg, r"\1", tree_newick, 1)      # Delete the outer parentheses from the analyzed clade
#        nodes.append([len(leaves), leaves, "[&&NHX:S=" + exp + ":D=N:B=speciesTree]"])
#
#    sortedNodes     = sorted(nodes, key=lambda x:x[0], reverse=True)
#
#    ## Add leaf as a clade
#    
#    largestClade = sortedNodes[0];
#    for leaf in largestClade[1]:
#        exp           = leaf
#        tempLeafClade = set([leaf])
#        sortedNodes.append([1, tempLeafClade, "[&&NHX:S=" + exp + "]"])
#
#    return sortedNodes


def collect_nodes_from_NHX(keyNode, treeFN):
    #print("### collect_nodes_from_NHX() ###")
    clades = []     # 2D array for clades
    expReg = r"\[.*?\]"
    cladeReg = r"\(([^\(\)]+)\)(.*?\[.*?\])"
    while re.search(cladeReg, treeFN):
        treeFN = treeFN.rstrip("\n")
        match = re.search(cladeReg, treeFN)                # Pick up the smallest and leftmost clade for the following analysis
        leavesString = match.group(1)
        exp = match.group(2)
        leavesString = re.sub(r"\[.*?\]", "", leavesString)   # Delete exp [...] of internal branches
        leavesString = re.sub(r":\d+\.\d+E-\d*", "", leavesString)   # Delete blanch lengths including E-
        leavesString = re.sub(r":\d+\.\d+", "", leavesString)       # Delete blanch lengths:  :0.013547
        leavesString = re.sub(r":-\d+\.\d+", "", leavesString) # Delete blanch lengths::-0.0
        #print("leavesString", leavesString)
        leaves = leavesString.split(",")
        treeFN = re.sub(cladeReg, r"\1", treeFN, count=1)                # Delete the outer parentheses from the analyzed clade
        #print("leaves", leaves)
        clades.append([len(leaves), leaves, exp])
    #print("### Exit point 1068 ###")
    #exit()

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


def identify_orthogroup(eachDirAddress, keyNode, treeNHX):
    print("### identify_orthogroup() ###")
    #print("keyNode", keyNode)
    #print("treeNHX", treeNHX)
    #f = open(eachDirAddress + "000_speciesTree_topLeft.txt")
    #speciesTree = "".join(list(f))
    #speciesTree = re.sub(r"[ \n]", "", speciesTree)
    #f.close()

    allGeneNodesSR = collect_nodes_from_NHX(keyNode, treeNHX)

    topHits = topHitPicker(eachDirAddress)
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
            #print("nameLine_cds_blastTopHit[1:]", nameLine_cds_blastTopHit[1:])
            if leaf == nameLine_cds_blastTopHit[1:]:
            #if re.search(nameLine_cds_blastTopHit[1:], leaf):
                criterion += 1
        #print("node[2]", node[2])
        if re.search(r"S=" + keyNode, node[2]):
            #print(" keyNode found")
            criterion += 1

        if criterion == 2:
           candidates_orthogroup.append(node)
           
        #print("   return")

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

    #print("orthogroupSR", orthogroupSR)
    #exit()
    return orthogroupSR


def leafCollectInOrderFrom_bothNHXnewick (tree):
    leaves = tree.split(",")
    leaves = [re.sub(r"^\(*([^:\(\)]+).*$", r"\1", leaf) for leaf in leaves]
    leaves = [re.sub(r"[\t\n]", "", leaf) for leaf in leaves]
    #for leaf in leaves:
    #    print(leaf)
    #exit()
    return leaves


def collect_speciesNames_in_orthogroup(SpeciesTree, keyNode, allNodes_speciesTree):
    #allNodes_speciesTree = collect_nodes_from_newick(SpeciesTree)
    leaves_species = leafCollectInOrderFrom_bothNHXnewick(SpeciesTree)
    orthoSpeciesGroup = ""
    speciesNames_in_orthogroup_FN = []
    for node in allNodes_speciesTree:
        if re.search(r"S=" + keyNode, node[2]):
            orthoSpeciesGroup = node
            break
    for leaf_species in leaves_species:
        if leaf_species in orthoSpeciesGroup[1]:
            speciesNames_in_orthogroup_FN.append(leaf_species)
    return speciesNames_in_orthogroup_FN


def identify_speciesNode(allNodes_speciesTree, name_node):
    for node in allNodes_speciesTree:
        if re.search(r"S=" + name_node, node[2]):
            return node


def identify_targetGeneNode(allNodes_speciesTree, allNodes_SR, name_speciesNode, separationType, queryGeneLeaf):
    #print("### identify_targetGeneNode() ###")
    #print("name_speciesNode", name_speciesNode)
    #print("separationType", separationType)
    #print("queryGeneLeaf", queryGeneLeaf)

    rec_targetSpeciesNode = identify_speciesNode(allNodes_speciesTree, name_speciesNode)
    #print("rec_targetSpeciesNode", rec_targetSpeciesNode)
    #exit()
    candidates_queryGeneNode = []

    flag_counting = 0
    for node in allNodes_SR:
        flag_counting += 1
        if flag_counting == 1:
            continue

        if separationType == "SisterGeneGroups" and flag_counting == 2:
            continue

        criterion = 0

        for leaf in node[1]:
            if leaf == queryGeneLeaf:
                criterion += 1

        keyWord = "S=" + name_speciesNode + r"[:\]]"
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
#            if re.search(r"S=" + name_focalGeneNode + ":", clade[2]):
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

def collect_AllchildSpeciesNodes_with_querySpecies(allspeciesNodesSR, focalNode_SR, name_querySpecies):
    childSpeciesNodes_with_querySpecies = []
    for speciesNode in allspeciesNodesSR:
        #print("speciesNode:", speciesNode[2])
        if name_querySpecies in speciesNode[1]:
            childSpeciesNodes_with_querySpecies.append(speciesNode)
    return childSpeciesNodes_with_querySpecies

def collect_childSpeciesNodes_with_querySpecies(allspeciesNodesSR, focalNode_SR, name_querySpecies):
    childSpeciesNodes_focalGroup = collect_childNodes(allspeciesNodesSR, focalNode_SR)
    childSpeciesNodes_with_querySpecies = []
    for speciesNode in childSpeciesNodes_focalGroup:
        if name_querySpecies in speciesNode[1]:
            childSpeciesNodes_with_querySpecies.append(speciesNode)
    return childSpeciesNodes_with_querySpecies


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
                match = re.search(r"S=([^:\]]+)[:\]]", branchLabel)
                hildBranchLavels.append(match.group(1))
        if re.search(r"S=" + nodeName + ":", node_speciesTree[2]):
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
#    match = re.search(r"S=([^:]+):", parentNode[2])
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


def identifiy_focalNode_speciesTree(allNodes_speciesTree, keyNode):
    #print("### identifiy_focalNode_speciesTree() ###")
    #print("keyNode", keyNode)
    focalNode_speciesTree = ""
    for node in allNodes_speciesTree:
        #print("node", node)
        if re.search(r"S=" + keyNode + ":", node[2]):
            focalNode_speciesTree = node
    #print("### Exit point")
    #exit()
    return focalNode_speciesTree


def count_duplications_for_speciesNodes(allGeneNodesSR, topHitName_1stQuery, childSpeciesNodes_focalGroup):
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
    #for childSpeciesNode_focalGroup in childSpeciesNodes_orthogorup:
    #    if not name_querySpecies in childSpeciesNode_focalGroup[1]:
    #        continue
    for childSpeciesNode_focalGroup in childSpeciesNodes_focalGroup:

        speciesNodeName = make_nodeName_from_nodeLavel_NHXstyle(childSpeciesNode_focalGroup[2])
        #print("speciesNodeName", speciesNodeName)
        count_dup = 0
        flag = 0
        for node_geneTree in geneNodes_including_querySequence:
            #print("node_geneTree[2]", node_geneTree[0], node_geneTree[2])
            if flag == 1:
                if re.search(r":D=Y[:\]]", node_geneTree[2]):
                    match = re.search(r"^([^\[]+)\[.*S=([^:]+):", node_geneTree[2])
                    nodeID_geneTree = match.group(1)
                    geneNodeName = match.group(2)
                    #print("nodeID_geneTree", nodeID_geneTree)
                    #print("geneNodeName", geneNodeName)
                    if geneNodeName == speciesNodeName and nodeID_geneTree.startswith("n"):
                        count_dup += 1
                        #print("count")

            #if flag == 1:
            #    if re.search(r"S=" + speciesNodeName + ":D=Y", node_geneTree[2]):
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
def count_blastHits(eachDirAddress, taxonSamplingList, infile):
    f = open(eachDirAddress + infile)
    lines = list(f)
    f.close()
    
    spNamePrefixes = []
    for spNameTMP in taxonSamplingList:
        match = re.search(r"^([^_]+_)([^_]+)$", spNameTMP)
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

    recs_cds_assigned_by_ID = readFasta_dict(eachDirAddress, "000_aaSeq_assigned_by_ID.txt")
    fs.write(">QuerySequence\n")
    lines_hit_query = make_lines_hit_query(eachDirAddress, recs_cds_assigned_by_ID)
    for line in lines_hit_query:
        #print("line QuerySequence:", line)
        fs.write(line)
    fs.write("\n")

    fs.close()


def make_bsvalue_from_nodeLavel_NHXstyle(nodeLavel_NHXstyle):
    #print("nodeLavel_NHXstyle", nodeLavel_NHXstyle)
    if nodeLavel_NHXstyle.startswith("["):
        return("leaf")
    match = re.search(r"B=([\d]+)", nodeLavel_NHXstyle)
    if match:
        return(match.group(1))
    else:
        return("r")

def make_duplicationStatus_from_nodeLavel_NHXstyle(nodeLavel_NHXstyle):
    #print("nodeLavel_NHXstyle", nodeLavel_NHXstyle)
    if nodeLavel_NHXstyle.startswith("["):
        return("leaf")
    if re.search(r":D=Y[:\]]", nodeLavel_NHXstyle):
        return("D=Y")
    if re.search(r":D=N[:\]]", nodeLavel_NHXstyle):
        return("D=N")
    else:
        print("Error in  make_duplicationStatus_from_nodeLavel_NHXstyle.")
        print("Check nodeLavel", nodeLavel_NHXstyle)
        exit()


def make_nodeName_from_nodeLavel_NHXstyle(nodeLavel_NHXstyle):
    #print("### make_nodeName_from_nodeLavel_NHXstyle() ###")
    #print("nodeLavel_NHXstyle", nodeLavel_NHXstyle)
    matchA = re.search(r"S=([^:]+)[:\]]", nodeLavel_NHXstyle)
    if matchA:
        return(matchA.group(1))
    else:
        return("No_node_lavel")


def check_sisterGroupName_included_in_ancestralSpeciesNodeNames(allNodes_speciesTree, allGeneNodesSR, targetSpeciesNode, sisterGeneGroup_SR):
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


def make_list_resLines_monophyletic(allNodes_speciesTree, allGeneNodesSR, nodes_speciesTree, topHitName_1stQuery):
    #print("### make_list_resLines_monophyletic() ####")
    list_resultLine = []
    #print("topHitName_1stQuery", topHitName_1stQuery)
    #exit()

    for node_speciesTree in nodes_speciesTree:

        name_speciesNode = make_nodeName_from_nodeLavel_NHXstyle(node_speciesTree[2])
        #print("name_speciesNode:", name_speciesNode)
        targetGeneNode = identify_targetGeneNode(allNodes_speciesTree, allGeneNodesSR, name_speciesNode, "MonophyleticGeneGroups", topHitName_1stQuery)
        #print("targetGeneNode:", targetGeneNode[0], targetGeneNode[2])
        whiteSpace = " " * (30 - len(name_speciesNode))
        if targetGeneNode[2] == "NoGeneNode":
            resultLine = name_speciesNode + "  " + whiteSpace + "NoGeneNode" + "  " + "NONE"
        else:
            resultLine = name_speciesNode + "  " + whiteSpace + make_bsvalue_from_nodeLavel_NHXstyle(targetGeneNode[2]) + "  " + make_duplicationStatus_from_nodeLavel_NHXstyle(targetGeneNode[2])
        #print("resultLine", resultLine)
        list_resultLine.append(resultLine)
    return list_resultLine


def make_list_resLines_sister(allNodes_speciesTree, allGeneNodesSR, nodes_speciesTree, topHitName_1stQuery):
    list_resultLine = []

    for node_speciesTree in nodes_speciesTree:

        name_speciesNode = make_nodeName_from_nodeLavel_NHXstyle(node_speciesTree[2])
        #print("name_speciesNode:", name_speciesNode)
        targetGeneNode = identify_targetGeneNode(allNodes_speciesTree, allGeneNodesSR, name_speciesNode, "SisterGeneGroups", topHitName_1stQuery)
        #print("targetGeneNode:", targetGeneNode[0], targetGeneNode[2])
        whiteSpace = " " * (30 - len(name_speciesNode))
        if targetGeneNode[2] == "NoGeneNode":
            resultLine = name_speciesNode + "  " + whiteSpace + "NONE   NoGeneNode"
        else:
            parentNode_queryGeneGroup, sisterGeneGroup = identify_sisterGeneNode(allGeneNodesSR, targetGeneNode)
            if sisterGeneGroup[2] == "NoSisterNode":
                resultLine = name_speciesNode + "  " + whiteSpace + "NONE   NoSisterNode"
            else:
                nodeName_sisterGeneGroup = check_sisterGroupName_included_in_ancestralSpeciesNodeNames(allNodes_speciesTree, allGeneNodesSR, node_speciesTree, sisterGeneGroup)
                resultLine = name_speciesNode + "  " + whiteSpace + make_bsvalue_from_nodeLavel_NHXstyle(parentNode_queryGeneGroup[2]) + "   " + nodeName_sisterGeneGroup
        list_resultLine.append(resultLine)
    return list_resultLine


#def add_makeSummary(ctx, outfile_summary2):
def add_makeSummary(
    eachDirAddress,
    treeSearchMethod,
    keyNode,
    allNodes_speciesTree,
    childSpeciesNodes_focalGroup,
    startTime,
    outfile_summary2,
):

    #print("#### add_makeSummary ####")
    fSum = open(eachDirAddress + outfile_summary2, "a")

    resDict_1stSummary = readRes_dict(eachDirAddress + "100_analysisSummary.txt")
    topHitName_1stQuery = resDict_1stSummary[">QuerySequence"][0]
    topHitName_1stQuery = re.sub(r" +.*", "", topHitName_1stQuery)
    #print("topHitName_1stQuery", topHitName_1stQuery)
    #exit()
    #topHitName_1stQuery = re.sub(r" .*$", "", topHitName_1stQuery)

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

    allGeneNodesSR_2ndTree = collect_nodes_from_NHX(keyNode, rearranged_2nd_gene_tree_NHX)

    fSum.write(">MonophyleticGeneGroups\n")
    list_resLines_mono = make_list_resLines_monophyletic(allNodes_speciesTree, allGeneNodesSR_2ndTree, childSpeciesNodes_focalGroup, topHitName_1stQuery)
    for line in list_resLines_mono:
        fSum.write(line + "\n")
    fSum.write("\n")

    #print("### >MonophyleticGeneGroups")
    ##for childSpeciesNode_focalGroup in childSpeciesNodes_orthogorup:
    ##    if not name_querySpecies in childSpeciesNode_focalGroup[1]:
    ##        continue
    #for childSpeciesNode in childSpeciesNodes_focalGroup:
    #
    #    name_speciesNode = make_nodeName_from_nodeLavel_NHXstyle(childSpeciesNode[2])
    #    print("name_speciesNode:", name_speciesNode)
    #    targetGeneNode = identify_targetGeneNode(allGeneNodesSR_2ndTree, name_speciesNode, "MonophyleticGeneGroups", topHitName_1stQuery)
    #    print("targetGeneNode:", targetGeneNode[0], targetGeneNode[2])
    #    whiteSpace = " " * (30 - len(name_speciesNode))
    #    if targetGeneNode[2] == "NoGeneNode":
    #        resultLine = name_speciesNode + "  " + whiteSpace + "NoGeneNode" + "  " + "NONE"
    #    else:
    #        resultLine = name_speciesNode + "  " + whiteSpace + make_bsvalue_from_nodeLavel_NHXstyle(targetGeneNode[2]) + "  " + make_duplicationStatus_from_nodeLavel_NHXstyle(targetGeneNode[2])
    #    print("resultLine", resultLine)
    #    fSum.write(resultLine + "\n")
    #fSum.write("\n")

    fSum.write(">SisterGeneGroups\n")
    list_resLines_sister = make_list_resLines_sister(allNodes_speciesTree, allGeneNodesSR_2ndTree, childSpeciesNodes_focalGroup, topHitName_1stQuery)
    for line in list_resLines_sister:
        fSum.write(line + "\n")
    fSum.write("\n")

    #fSum.write(">SisterGeneGroups\n")
    #print("### >SisterGeneGroups")
    ##for targetSpeciesNode in speciesNodes_including_querySpecies:
    ##for childSpeciesNode_focalGroup in childSpeciesNodes_orthogorup:
    ##    if not name_querySpecies in childSpeciesNode_focalGroup[1]:
    ##        continue
    #for ChildSpeciesNode in childSpeciesNodes_focalGroup:
    #
    #    name_speciesNode = make_nodeName_from_nodeLavel_NHXstyle(ChildSpeciesNode[2])
    #    print("name_speciesNode:", name_speciesNode)
    #    targetGeneNode = identify_targetGeneNode(allGeneNodesSR_2ndTree, name_speciesNode, "SisterGeneGroups", topHitName_1stQuery)
    #    print("targetGeneNode:", targetGeneNode[0], targetGeneNode[2])
    #    whiteSpace = " " * (30 - len(name_speciesNode))
    #    if targetGeneNode[2] == "NoGeneNode":
    #        resultLine = name_speciesNode + "  " + whiteSpace + "NONE   NoGeneNode"
    #    else:
    #        parentNode_queryGeneGroup, sisterGeneGroup = identify_sisterGeneNode(allGeneNodesSR_2ndTree, targetGeneNode)
    #        if sisterGeneGroup[2] == "NoSisterNode":
    #            resultLine = name_speciesNode + "  " + whiteSpace + "NONE   NoSisterNode"
    #        else:
    #            nodeName_sisterGeneGroup = check_sisterGroupName_included_in_ancestralSpeciesNodeNames(allGeneNodesSR_2ndTree, ChildSpeciesNode, sisterGeneGroup)
    #            resultLine = name_speciesNode + "  " + whiteSpace + make_bsvalue_from_nodeLavel_NHXstyle(parentNode_queryGeneGroup[2]) + "   " + nodeName_sisterGeneGroup
    #    print("resultLine", resultLine)
    #    fSum.write(resultLine + "\n")
    #fSum.write("\n")
    ####

    #fSum.write(">BootstrapValue_sisterGeneGroup\n")
    #fSum.write(make_bsvalue_from_nodeLavel_NHXstyle(sisterGeneGroup[2]) + "\n\n")
    #fSum.write(">Members_sisterGeneGroup\n")
    #for member in sisterGeneGroup[1]:
    #    fSum.write(member + "\n")
    #fSum.write("\n")

    fSum.write(">Number_of_duplicatedNode\n")
    recs_duplications_for_speciesNodes = count_duplications_for_speciesNodes(allGeneNodesSR_2ndTree, topHitName_1stQuery, childSpeciesNodes_focalGroup)
    for nodeName, numDup in recs_duplications_for_speciesNodes.items():
        whiteSpace =  " " * (30 - len(nodeName)) 
        fSum.write(nodeName + "  " + whiteSpace + str(numDup) + "\n")
    fSum.write("\n")

    fSum.write(">AnalysisTime\n")
    elapsed_time = round((time.time() - startTime),1)
    fSum.write (str(elapsed_time) + " seconds")
    fSum.write("\n")

    fSum.close()


def makeSummary(
    SpeciesTree,
    taxonSamplingList,
    mode,
    dataset,
    keyNode,
    startTime,
    eachDirAddress,
    eachDirAddress_e1stre,
    BSthreshold,
    BSthreshold_4_ReE1st,
    num_rootSequences,
    allNodes_speciesTree,
    childSpeciesNodes_AllGroup,
    *,
    aligned_site_rate=None,
    outfile_summary="100_analysisSummary.txt",
    fatal_error_msg=None,   # ← 追加
):
    print("#### makeSummary ####")

    # 追加：致命的エラー時は、最低限の Summary を作って返す
    if fatal_error_msg:
        fs = open(eachDirAddress + outfile_summary, "w")
        fs.write("################ Results: 1st analysis ################\n\n")
        fs.write(">BS_of_orthogroupBasalNode\n")
        fs.write(fatal_error_msg + "\n\n")

        fs.write("\n################ Settings ################\n\n")
        fs.write(">Mode\n")
        fs.write(mode + "\n\n")
        fs.write(">Dataset\n")
        fs.write(dataset + "\n\n")
        fs.write(">ShortSequence_threshold\n")
        fs.write(str(aligned_site_rate) + "\n\n")
        fs.write(">SpeciesTree\n")
        fs.write(SpeciesTree + "\n\n")

        fs.write(">AnalysisTime\n")
        elapsed_time = round((time.time() - startTime), 1)
        fs.write(str(elapsed_time) + " seconds\n")
        fs.close()
        return

    recs_cds_assigned_by_ID = readFasta_dict(eachDirAddress, "000_aaSeq_assigned_by_ID.txt")
    #cDNAfn = ""
    #cDNAfn = readFasta_dict(eachDirAddress, "000_aaSeq_assigned_by_ID.txt")
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
    lines_query = make_lines_hit_query(eachDirAddress, recs_cds_assigned_by_ID)
    for line in lines_query:
        #print("line:", line)
        fs.write(line)
    fs.write("\n")

    topHitName_1stQuery, dummy1 = lines_query[0].split("<=")
    topHitName_1stQuery = re.sub(" *$", "", topHitName_1stQuery)
    #print("topHitName_1stQuery:", topHitName_1stQuery, "|")

    fs.write(">Number_of_blastHits\n")
    rec_blastHits = count_blastHits(eachDirAddress, taxonSamplingList, infile = "010_blastRes.txt")
    #print("rec_blastHits", rec_blastHits)
    #exit()
    rec_blastHits = whiteSpaceAdd(rec_blastHits)
    for name, num in rec_blastHits.items():
        fs.write(name + str(num) + "\n")
    fs.write("\n")

    #speciesTree = ""
    #f1stTree = open(eachDirAddress + "000_speciesTree_topLeft.txt")
    #speciesTree = "".join(list(f1stTree))
    #speciesTree = re.sub(r"[ \n]", "", speciesTree)
    #f1stTree.close()

    rearranged_1st_gene_tree_NHX = ""
    fname = os.path.join(eachDirAddress, "080_trimedAAPhy.txt")
    if not os.path.isfile(fname):
        print("No file:", fname)
        outgroup1 = "No_file"
    else:
        outgroup1 = outGroupSelect(eachDirAddress, "080_trimedAAPhy.txt")
    path = os.path.join(eachDirAddress, "085_NJBS1st.txt.rearrange.0")
    if os.path.exists(path):
        fs.write(">Rooting_4_1stTree\n")
        fs.write(f"{outgroup1}\n")
        fs.write("\n")
        with open(path) as f1stTree:
            rearranged_1st_gene_tree_NHX = f1stTree.readline().rstrip("\n")

    if not fMafOut:
        fs.write(">BS_of_orthogroupBasalNode\n")
        fs.write("No mafft out.\n")
        fs.write("\n")
    elif not rearranged_1st_gene_tree_NHX:
        fs.write(">BS_of_orthogroupBasalNode\n")
        fs.write("1st tree not estimated.\n")
        fs.write("\n")
    else:
        orthogroup = identify_orthogroup(eachDirAddress, keyNode, rearranged_1st_gene_tree_NHX)
        if orthogroup[0] == 0:
            fs.write(">BS_of_orthogroupBasalNode\n")
            fs.write(orthogroup[2] + "\n")
            fs.write("\n")
        else:
            fs.write(">BS_of_orthogroupBasalNode\n")
            if len(orthogroup[1]) < 4:
                fs.write("Less than 4 orthogroup members.\n")

            elif mode == "E1st":
                fs.write("1st tree estimated by mode E1st.\n")
            else:
                fs.write(make_bsvalue_from_nodeLavel_NHXstyle(orthogroup[2]) + "\n\n")
            fs.write("\n")

            if mode == "E":
                fs.write(">Orthogroup\n")
                sorted_members_focalGeneNode = sorted(orthogroup[1])
                for leaf in sorted_members_focalGeneNode:
                    fs.write(leaf + "\n")
                fs.write("\n")
        
                fs.write(">GeneNumber_of_orthogroup\n")
                speciesNames_in_orthogroup = collect_speciesNames_in_orthogroup(SpeciesTree, keyNode, allNodes_speciesTree)
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

            fs.write(">Rooting_4_2ndTree\n")
            rootGeneLeaves = selectRootSp4secondTreeSearch(eachDirAddress, keyNode, num_rootSequences)
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


        ### START: Identifying monophyletic/sister gene nodes
        allGeneNodesSR_1stTree = collect_nodes_from_NHX(keyNode, rearranged_1st_gene_tree_NHX)
    
        fs.write(">MonophyleticGeneGroups_1stTree\n")
        list_resLines_mono = make_list_resLines_monophyletic(allNodes_speciesTree, allGeneNodesSR_1stTree, childSpeciesNodes_AllGroup, topHitName_1stQuery)
        for line in list_resLines_mono:
            fs.write(line + "\n")
        fs.write("\n")
    
        fs.write(">SisterGeneGroups_1stTree\n")
        list_resLines_sister = make_list_resLines_sister(allNodes_speciesTree, allGeneNodesSR_1stTree, childSpeciesNodes_AllGroup, topHitName_1stQuery)
        for line in list_resLines_sister:
            fs.write(line + "\n")
        fs.write("\n")
        ### END: Identifying monophyletic/sister gene nodes


        if rec044_unambSiteRate:
            rec044_unambSiteRate = whiteSpaceAdd(rec044_unambSiteRate)
            fs.write(">Aligned-ShortSequence_threshold evaluation\n")
            for name, siteRate in rec044_unambSiteRate.items():
                nameRate = name + "  " + str(siteRate)
                if float(siteRate) > float(aligned_site_rate):
                    fs.write(nameRate + "\n")
                else:
                    fs.write("== Removed ==> " + nameRate + "\n")
            fs.write("\n")


    fs.write("\n################ Settings ################\n\n")

    fs.write(">Mode\n")
    fs.write(mode + "\n")
    fs.write("\n")

    fs.write(">NumberAssigned_querySequence\n")
    for name, seq in recs_cds_assigned_by_ID.items():
        name = re.sub(r"[\n\r]", "", name)
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

    fs.write(">ShortSequence_threshold\n" + str(aligned_site_rate) + "\n\n")

    fs.write(">Dataset\n" + dataset +  "\n")
    fs.write("\n")

    fs.write(">Outdir\n" + dataset +  "\n")
    fs.write("\n")

    fs.write(">Outdir_4_ReE1st\n" + dataset +  "\n")
    fs.write("\n")

    fs.write(">Rearrangement_BS_value_threshold\n")
    fs.write(str(BSthreshold) + "\n")
    fs.write("\n")

    fs.write(">BSthreshold_4_ReE1st\n")
    fs.write(str(BSthreshold_4_ReE1st) + "\n")
    fs.write("\n")


    fs.write(">TaxonSampling_color\n")
    for spNameTMP in taxonSamplingList:
        fs.write(spNameTMP + "\n")
    fs.write("\n")

    fs.write(">AnalysisTime\n")
    elapsed_time = round((time.time() - startTime),1)
    fs.write (str(elapsed_time) + " seconds")
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

def select_blastHitUnique(eachDirAddress, blastResFileFN):
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
            nameLine_tmp = re.sub(r"^> ", ">", nameLine_tmp)
            match = re.search(r"^ (Identities = [^,]+),", line)
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
    #nameLine_FN = re.sub(r"\+\+\-", "-", nameLine_FN)
    nameLine_FN = re.sub(r"\+\+", "+", nameLine_FN)
    nameLine_FN = re.sub(r"-+", "-", nameLine_FN)
    return nameLine_FN

def shorten_nameLine(nameLineTMP):
    #print("### shorten_nameLine ###")
    nameLine = re.sub(r"(>.{60}).*", r"\1", nameLineTMP)
    #print("nameLine", nameLine)
    return nameLine

def hitRecPicker(dbAddress, dbLines, eachDirAddress):
    #print("#### hitRecPicker ####")
    blastResOut = open(eachDirAddress + "010_blastRes.txt", "w")
    AAout = open(eachDirAddress + "030_retrievedAAfas.txt",  "w")
    CDNAout = ""
    CDNAout = open(eachDirAddress + "030_retrievedDNAfas.txt", "w")

    num_BlastpHits = 0
    for dbline in dbLines:

        recdbAA = readFasta_dict(dbAddress, dbline[1])
        recdbDNA = ""
        recdbDNA = readFasta_dict(dbAddress, dbline[2])

        blastResFile   = "005_vs" + dbline[0][:-1] + ".txt"
        #print("blastResFile:", blastResFile)
        recs_nameLine_identity_unique = select_blastHitUnique(eachDirAddress, blastResFile)
        
        for nameLine_blasthit in reversed(recs_nameLine_identity_unique.keys()):
            if not nameLine_blasthit.startswith(">"):
                continue
            
            num_BlastpHits += 1

            nameLine_blasthit = nameLine_blasthit.rstrip("\n")
            nameLine_blasthitTMP = re.sub(r">", ">" + dbline[0], nameLine_blasthit)
            blastResOut.write(nameLine_blasthitTMP + "\n")

            counter_DBNLINEnum = 0
            for name_recdbAA, seq_recdbAA in recdbAA.items():
                if name_recdbAA == nameLine_blasthit:
                    break
                counter_DBNLINEnum += 1
            #match = re.search(r"DBNLINE\|(\d+)\|", nameLine_blasthit)
            #DBNLINEnum = match.group(1)
            DBNLINEnum = counter_DBNLINEnum
            nameLine = list(recdbAA.keys())[int(DBNLINEnum)]
            nameLine = change_prohibitedExpression_in_nameLine(nameLine)
            nameLine = re.sub(r">", ">" + dbline[0][:-1] + "_", nameLine)
            #print("nameLine",nameLine)
            
            #print("nameLine1", nameLine)
            nameLine = re.sub(" .*$", "", nameLine)
            #print("nameLine2", nameLine)
            nameLine = shorten_nameLine(nameLine)
            #print("nameLine3", nameLine)
            #print("")

            AAout.write(nameLine + "\n")
            AAout.write(list(recdbAA.values())[int(DBNLINEnum)] + "\n")
            #CDNAout.write(list(recdbDNA.keys())[int(DBNLINEnum)]   + "\n")
            CDNAout.write(nameLine + "\n") # modified 20241118
            CDNAout.write(list(recdbDNA.values())[int(DBNLINEnum)] + "\n")
            
            #print()
    blastResOut.close()

    if num_BlastpHits < 4:
        result = "Less than 4 blast hits."
        error_makeSummary(result)
        #print (result)
        if Switch_deleteIntermediateFiles == "L":
            error_resHtmlMaker(ctx.eachDirAddress, ctx.keyNode, ctx.queryID, result)
        if Switch_deleteIntermediateFiles == "D":
            deleteFiles(eachDirAddress)
        exit()

    AAout.close()
    CDNAout.close()


def make_spNames_IndbLines():
    spNames_IndbLines = []
    for dbLine in dbLines:
        spNamePrefix_dbLine = dbLine[0]
        spNamePrefix_dbLine = re.sub(r"_$", "", spNamePrefix_dbLine)
        spNames_IndbLines.append(spNamePrefix_dbLine)
    return spNames_IndbLines


def blastpSearch(dbLines, dbAddress, toolAddress, blastEvalue, Number_of_hits_to_report_per_genome, eachDirAddress):
    DirAafileName = eachDirAddress + "000_translated_cds_assigned_by_ID.txt"
    for dbline in dbLines:
        dbAAfile = dbAddress + dbline[1]
        outFile = eachDirAddress + "005_vs" + dbline[0][:-1] + ".txt"
        #comLine = "tools/blastp -query {0} -evalue {1} -num_alignments {2}  -num_descriptions {3}  -db {4} -out {5}"\
        #          .format(DirAafileName, blastEvalue, Number_of_hits_to_report_per_genome, Number_of_hits_to_report_per_genome, dbAAfile, outFile)
        comLine = toolAddress + "blastp -query {0} -evalue {1} -num_alignments {2}  -num_descriptions {3}  -db {4} -out {5}"\
                  .format(DirAafileName, blastEvalue, Number_of_hits_to_report_per_genome, Number_of_hits_to_report_per_genome, dbAAfile, outFile)
        #print("comLine", comLine)
        subprocess.call(comLine, shell=True)


def make_lines_hit_query(eachDirAddress, recs_cds_assigned_by_ID):
    topHits = topHitPicker(eachDirAddress)

    topHitValue0s = [x[0] for x in topHits.values()]
    length_longestName = len(max(topHitValue0s, key = len))

    lines_hit_query_FN = []
    for i in range(len(topHits)):
        uploadedSeqName = list(recs_cds_assigned_by_ID.keys())[i][1:]
        uploadedSeqName = re.sub(r"[\n\r]", "", uploadedSeqName)
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


def topHitPicker(eachDirAddress):
    topHits = OrderedDict()
    recs_cds_assigned_by_ID = readFasta_dict(eachDirAddress, "000_aaSeq_assigned_by_ID.txt")
    for nameLine_cds_assigned_by_ID in recs_cds_assigned_by_ID.keys():
        match = re.search(r">([^_]+)_", nameLine_cds_assigned_by_ID)
        querySpecies_SR = match.group(1)
        #recs_blastnRes = read_blastnRes2(eachDirAddress, "005_vs" +querySpecies_SR + ".txt")
        recs_blastnRes = select_blastHitUnique(eachDirAddress, "005_vs" + querySpecies_SR + ".txt")
        blastTopHitNameLine = ""
        blastTopHitIdentity = ""
        for nameLine_blasthit, identity_blasthit in recs_blastnRes.items():
            blastTopHitNameLine = re.sub(r" .*$", "", nameLine_blasthit)
            #print("blastTopHitNameLine", blastTopHitNameLine)
            blastTopHitNameLine = change_prohibitedExpression_in_nameLine(blastTopHitNameLine)
            #print("blastTopHitNameLine", blastTopHitNameLine)
            print("")
            blastTopHitNameLine = re.sub(r">", ">" + querySpecies_SR + "_", blastTopHitNameLine)
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
        if re.match(r"\w", querySeq[i]) and re.match(r"\w", otherSeq[i]):
            count += 1
        elif re.match(r"-", querySeq[i]) and re.match(r"\w", otherSeq[i]):
            count += 1
        else:
            continue
    return round(float(count)/len(querySeq),3)

def compare_query2other_aliSiteRate_nogap(querySeq, otherSeq):
    count = 0
    for i in range(len(querySeq)):
        if re.match(r"-", querySeq[i]):
            continue
        if re.match(r"\w", otherSeq[i]):
            count += 1
        else:
            continue
    querySeq_noGap = re.sub(r"-", "", querySeq)
    #print("querySeq", querySeq)
    #print("querySeq_noGap", querySeq_noGap)
    #exit()
    return round(float(count)/len(querySeq_noGap),3)


def calculate_nonGapSiteRate(eachDirAddress, trimledAAfile):
    #print("#### calculate_nonGapSiteRate ####")
    ## trimledAAfile 042_AA.fas.trm
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
        protID = re.sub(r" .*$","", otherName)
        #aaRate = compare_query2other_aliSiteRate(querySeq, otherSeq)
        aaRate = compare_query2other_aliSiteRate_nogap(querySeq, otherSeq)
        #print("protID", protID)
        #print("aaRate", aaRate)
        #print("querySeq", querySeq)
        #print("otherSeq", otherSeq)
        recs_nonGapSiteRateFN[protID] = aaRate
    #exit()
    return recs_nonGapSiteRateFN


def delete_sequences_with_alignedSiteRate(eachDirAddress, aligned_site_rate, outFile_AA, outFile_DNA, outFile_siteRate):
    print("### delete_sequences_with_alignedSiteRate() ###")
    #print("outFile_AA", outFile_AA)
    #print("outFile_DNA", outFile_DNA)
    #print("outFile_siteRate", outFile_siteRate)

    dic_nonGapSiteRate = calculate_nonGapSiteRate(eachDirAddress, "042_AA.fas.trm")
    #for name, rate in dic_nonGapSiteRate.items():
    #    print("name1", name)
    #    print("rate1", rate)
    #exit()
    #print("### delete_sequences_with_alignedSiteRate 2 ###")

    file_overRateAA = open(eachDirAddress + outFile_AA, "w")
    file_overRateDNA = open(eachDirAddress + outFile_DNA, "w")
    file_unambSiteRatefile = open(eachDirAddress + outFile_siteRate, "w")

    dic_AA = readFasta_dict(eachDirAddress, "030_retrievedAAfas.txt")
    dic_DNA = ""
    dic_DNA = readFasta_dict(eachDirAddress, "030_retrievedDNAfas.txt")

    #print("len(dic_nonGapSiteRate)", len(dic_nonGapSiteRate))
    #print("len(dic_AA)",len(dic_AA))
    #print("len(dic_DNA)",len(dic_DNA))
    #exit()


    #if len(dic_nonGapSiteRate) == len(dic_AA):
    #    for i in range (len(dic_nonGapSiteRate)):
    #        name = list(dic_nonGapSiteRate.keys())[i]
    #        rate = list(dic_nonGapSiteRate.values())[i]
    #        #print("name1", name)
    #        #print("rate1", rate)
    #        #print("aligned_site_rate", float(aligned_site_rate))
    #        if rate > float(aligned_site_rate):
    #            print(rate, aligned_site_rate)
    #            #print("(list(dic_AA.keys())[i]", i, list(dic_AA.keys())[i]) # modified 20241118
    #            file_overRateAA.write(list(dic_AA.keys())[i] + "\n")
    #            #print("list(dic_AA.keys())[i]", list(dic_AA.keys())[i])
    #            file_overRateAA.write(list(dic_AA.values())[i] + "\n")
    # 
    #            #print("(list(dic_DNA.keys())[i]", i, list(dic_DNA.keys())[i])
    #            file_overRateDNA.write(list(dic_DNA.keys())[i] + "\n")
    #            #print("list(dic_DNA.keys())[i]", list(dic_DNA.keys())[i])
    #            file_overRateDNA.write(list(dic_DNA.values())[i] + "\n")
    #
    #        file_unambSiteRatefile.write(name + "\n" + str(rate) + "\n")
    #        #print("")
    #else:
    #    ### 042_AA.fas.trm  : Only gap sequence is moved by TRIMAL1.4
    #    for name_siterate, rate_siterate in dic_nonGapSiteRate.items():
    #        #print("name_siterate", name_siterate)
    #        match_nsr = re.search(r"^[^_]+_([^_]+)_", name_siterate)
    #        name_siterate_geneID = match_nsr.group(1)
    #        #print("name_siterate_geneID", name_siterate_geneID)
    #        #print("rate_siterate", rate_siterate)
    #        if float(rate_siterate) > float(aligned_site_rate):
    #            
    #            for name_AA, seq_AA in dic_AA.items():
    #                #print("name_AA", name_AA)
    #                if re.search(r"_" + name_siterate_geneID + "_", name_AA):
    #                    file_overRateAA.write(name_AA + "\n")
    #                    file_overRateAA.write(seq_AA + "\n")
    #                    break
    #
    #            for name_DNA, seq_DNA in dic_DNA.items():
    #                #print("name_DNA", name_DNA)
    #                #if re.search(r">" + name_siterate_geneID + " ", name_DNA):
    #                if re.search(r"_" + name_siterate_geneID + "_", name_DNA):
    #                    file_overRateDNA.write(name_DNA + "\n")
    #                    file_overRateDNA.write(seq_DNA + "\n")
    #                    break
    #
    #        file_unambSiteRatefile.write(name_siterate + "\n" + str(rate_siterate) + "\n")
    #        #print("")

    for name_siterate, rate_siterate in dic_nonGapSiteRate.items():
        #print("name_siterate")
        #print("|", name_siterate, "|")
        #match_nsr = re.search(r"^[^_]+_([^_]+)_", name_siterate)
        #name_siterate_geneID = match_nsr.group(1)
        #print("name_siterate_geneID", name_siterate_geneID)
        #print("rate_siterate", rate_siterate)
        if float(rate_siterate) > float(aligned_site_rate):
            
            #for name_AA, seq_AA in dic_AA.items():
            #    print("name_AA")
            #    print("|", name_AA, "|")
            #    #if re.search(r"_" + name_siterate_geneID + "_", name_AA):
            #    #    #print("name_AA", name_AA)
            #    #    file_overRateAA.write(name_AA + "\n")
            #    #    file_overRateAA.write(seq_AA + "\n")
            #    #     break
            #    if name_siterate.rstrip("\n") == name_AA.rstrip("\n"):
            #        print("same")
            #    else:
            #        print("diff")
            #    #print("")
            file_overRateAA.write(name_siterate + "\n")
            file_overRateAA.write(dic_AA[name_siterate] + "\n")

            #for name_DNA, seq_DNA in dic_DNA.items():
            #    #print("name_DNA", name_DNA)
            #    #if re.search(r">" + name_siterate_geneID + " ", name_DNA):
            #    if re.search(r"_" + name_siterate_geneID + "_", name_DNA):
            #        file_overRateDNA.write(name_DNA + "\n")
            #        file_overRateDNA.write(seq_DNA + "\n")
            #        break
            file_overRateDNA.write(name_siterate + "\n")
            file_overRateDNA.write(dic_DNA[name_siterate] + "\n")

        file_unambSiteRatefile.write(name_siterate + "\n" + str(rate_siterate) + "\n")
        #print("")

    file_overRateAA.close()
    file_overRateDNA.close()
    file_unambSiteRatefile.close()
    #exit()

def compare_numSeqs(eachDirAddress, file_mafOutAA, file_overRateAA):
    #print("### compare_numSeqs() ###")
    #print("eachDirAddress", eachDirAddress)
    #print("file_mafOutAA", file_mafOutAA)
    #print("file_overRateAA", file_overRateAA)
    dic_mafOutAA = readFasta_dict(eachDirAddress, file_mafOutAA)
    dic_overRateAA = readFasta_dict(eachDirAddress, file_overRateAA)
    #for name, seq in dic_overRateAA.items():
    #    print("name", name)
    #    print("seq", seq)
    if len(dic_mafOutAA) == len(dic_overRateAA):
        return "Equal"
    else:
        return "NotEqual"
    exit()
    

### delete_sequences_with_alignedSiteRate End
###############################################################


###############################################################
### alignmentFile_PhyAnal Start
def orderedDict2FasFile(eachDirAddress, recs, outfile):
    out = open(eachDirAddress + outfile, "w")
    for name,value in recs.items():
        out.write(name + "\n")
        out.write(value + "\n")
    out.close()


def orderedDict2phyFile(eachDirAddress, recs, outfile):
    secLength = len(sorted(recs.values())[0])
    spSeqSizeLine = str(len(recs)) + " " + str(secLength)

    recs = whiteSpaceAdd(recs)
    out = open(eachDirAddress + outfile, "w")
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
        if re.search(r"    <span class=sel>Selected Sequences",line):
            continue
        if line.startswith("    <span class=sel>"):
            line = line.rstrip("\n")
            match = re.search(r"<span class=sel>([^<]+)</span> +([^ ].*)$", line)
            name     = match.group(1)
            sequence = match.group(2)
            if not name in recs_trimalHTMLout.keys():
                recs_trimalHTMLout[name] = ""
            recs_trimalHTMLout[name] += sequence

    return recs_trimalHTMLout

def read_TrimalHTMLout_dict_v141(eachDirAddress, trimAAresult):
    f = open(eachDirAddress + trimAAresult)
    lines = list(f)
    f.close()

    trimalMarkedSites = ""
    for line in lines:
        if line.startswith("    Selected Cols:"):
            line = line.rstrip("\n")
            sequence = re.sub(r" +Selected Cols: +", "",line)
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
    trimalMarkedSites = re.sub(r"<span class=sel>.</span>", "#", trimalMarkedSites)
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


def trimaledv141_FileMakerDNA(eachDirAddress, fastaFile, trimAAresult, outfile):
    recs = readFasta_dict(eachDirAddress, fastaFile)
    #print("fastaFile", fastaFile)
    #exit()
    
    trimalMarkedSitesTMP = read_TrimalHTMLout_dict_v141(eachDirAddress, trimAAresult)
    #print("trimalMarkedSitesTMP", trimalMarkedSitesTMP)
    #exit()
    trimalMarkedSites = re.sub(r"<span class=nsel> </span>", "-", trimalMarkedSitesTMP)
    trimalMarkedSites = re.sub(r"<span class=sel> </span>", "#", trimalMarkedSites)
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
    orderedDict2phyFile(eachDirAddress, recsTrimed, outfile = outfile)


### alignmentFile_PhyAnal End
###############################################################


###############################################################
### alignmentFile_HTML
def outGroupSelect(eachDirAddress, phyFileName):
    recSeqFN = readPhy_dict(eachDirAddress, phyFileName)
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
    recsFN2 = OrderedDict()
    for name, seq in recsFN.items():
        if re.search(r" \d+ bp$", name):
            name = re.sub(r" \d+ bp$", "", name)
        recsFN2[name] = seq
    return recsFN2


def nameChange_whiteLaterDelete(recsFN):
    recsFN2        = OrderedDict()
    for name, seq in recsFN.items():
        name = re.sub(r" .*$", "", name)
        recsFN2[name] = seq
    return recsFN2


def gapDelete(recsFN):
    recsFN2        = OrderedDict()
    for name, seq in recsFN.items():
        seq = re.sub(r"-", "", seq)
        recsFN2[name] = seq
    return recsFN2


def whiteSpaceAdd(recsFN1):
    longestName = max(recsFN1.keys(), key = len)
    longestName = re.sub(r"<[^>]+>", "", longestName)
    longestNameLen = len(longestName)
    recsFN2        = OrderedDict()
    for name,value in recsFN1.items():
        name = re.sub(r"^>", "", name)
        nameTMP = re.sub(r"<[^>]+>", "", name)
        nameWhiteSpace = name + " " * (longestNameLen - len(nameTMP) + 2)
        recsFN2[nameWhiteSpace] = value
    return recsFN2


def reorderDeleteGap_Fas2FasByTree(fastaFileName, tree4leafOrder, outPhyFileName):
    recs = readFasta_dict(eachDirAddress,fastaFileName) 
    recs = nameChange_whiteLaterDelete(recs)
    recs = gapDelete(recs)
    recs = reorderSeqByTree(recs, tree4leafOrder)
    orderedDict2FasFile(recs, outfile = outPhyFileName)


def phy2fastmePhy(eachDirAddress, phyFileName, outFastmePhyFileName):
    recsFN = readPhy_dict(eachDirAddress, phyFileName)
    secLength = len(sorted(recsFN.values())[0])
    spSeqSizeLine = str(len(recsFN)) + " " + str(secLength)

    recsFN = whiteSpaceAdd(recsFN)
    out = open(eachDirAddress + outFastmePhyFileName, "w")
    out.write(spSeqSizeLine + "\n")
    for name,value in reversed(recsFN.items()):
        out.write(name + value + "\n")
    out.close()


def fas2phy(eachDirAddress, fastaFileName, outPhyFileName):
    recs = readFasta_dict(eachDirAddress, fastaFileName)
    recs = delete_nameSpaceSeqBp(recs)
    orderedDict2phyFile(eachDirAddress, recs, outfile = outPhyFileName)

def phy2fas(eachDirAddress, infile_phy, outfile_fas):
    recs = readPhy_dict(eachDirAddress, infile_phy)  # PHYLIPファイルを読み込む
    orderedDict2FasFile(eachDirAddress, recs, outfile_fas)

def reorderFas2PhyByTree(fastaFileName, tree4leafOrder, outPhyFileName):
    print("### reorderFas2PhyByTree() ###")
    recs = readFasta_dict(eachDirAddress,fastaFileName) 
    recs = delete_nameSpaceSeqBp(recs)
    for name, seq in recs.items():
        print(name)
        print("")
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


def phyCodonToBlock(eachDirAddress, phyFileName, blockNum, outfile):
    recs = readPhy_dict(eachDirAddress, phyFileName)
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


def selectRootSp4secondTreeSearch(eachDirAddress, keyNode, num_rootSequences):
    f1stTree = open(eachDirAddress + "085_NJBS1st.txt.rearrange.0")
    treeNHX = list(f1stTree)[0]
    f1stTree.close()
    allGeneNodes = collect_nodes_from_NHX(keyNode, treeNHX)
    orthogroup = identify_orthogroup(eachDirAddress, keyNode, treeNHX)
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
            topHits = topHitPicker(eachDirAddress)
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


def make_2ndanalysis_seqFile(eachDirAddress, rootLeaves_SR, seqfile, resDict_1st, outfile):
    recs_nucl = readFasta_dict(eachDirAddress, seqfile)
    leaves_add_2_focalClade = []
    for rootLeaf in rootLeaves_SR:
       if not rootLeaf in resDict_1st[">Orthogroup"]:
            leaves_add_2_focalClade.append(rootLeaf)
    for leaf in resDict_1st[">Orthogroup"]:
        leaves_add_2_focalClade.append(leaf)

    out = open(eachDirAddress + "/" + outfile, "w")
    for name in leaves_add_2_focalClade:
        out.write(">" + name + "\n")
        seq = re.sub(r"-", "", recs_nucl[">" + name])
        out.write(seq + "\n")
    out.close()


def cDNAfas2noGapAAFasFile(eachDirAddress, cDNAfasFileName, outfile):
    recsFN = readFasta_dict(eachDirAddress, cDNAfasFileName)
    fa = open(eachDirAddress + "/" + outfile, "w")
    for name, seq in recsFN.items():
        fa.write(name + "\n")
        fa.write(translation(seq) + "\n")
    fa.close()


def make_raxmlPartitionFile(eachDirAddress, outPartFile):
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


def moveRAxMLfiles(eachDirAddress, outfile):
    line1 = "cp " + eachDirAddress + "RAxML_bipartitions.txt " + eachDirAddress + outfile
    #print(line1)
    subprocess.call(line1, shell=True)
    line2 = "rm " + eachDirAddress + "RAxML*"
    subprocess.call(line2, shell=True)

### alignmentFile_HTML End
###############################################################


###############################################################
### result_HTML Start
def error_resHtmlMaker(eachDirAddress, keyNode, queryID, resultFN):
    topHits = topHitPicker(eachDirAddress)
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


def make_resHtml2(
    queryID,
    mode,
    BSthreshold_4_ReE1st,
    eachDirAddress_FN,
    resHTMLlinesFN
):
    #print("### make_resHtml2() ###")
    #print("eachDirAddress_FN", eachDirAddress_FN)
    #exit()
    address_file_summary = os.path.join(eachDirAddress_FN, "100_analysisSummary.txt")
    resDict_1st = readRes_dict(address_file_summary)
    topHitName_1stQuery = resDict_1st[">QuerySequence"][0]
    topHitName_1stQuery = re.sub(r" +.*", "", topHitName_1stQuery)
    firstQueryTMP = shorten_nameLine(topHitName_1stQuery)
    resHTMLlines1 = re.sub('FIRSTQUERY', topHitName_1stQuery, resHTMLlinesFN)

    #f1stTree = open(eachDirAddress_FN + "085_NJBS1st.txt.rearrange.0")
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

    resHTMLlines1 = re.sub('EACHDIRADDRESS_', eachDirAddress_FN, resHTMLlines1)

    #print("mode", mode)
    #print("BSthreshold_4_ReE1st", BSthreshold_4_ReE1st)
    if mode in ("E1st", "D1st", "ReE1st"):
        resHTMLlines1 = re.sub('<td colspan="2"><b>2nd tree:</b> Speciation/duplication events in the query sequence lineage </td>', '<td colspan="2"><b>2nd tree:</b> Not estimated </td>', resHTMLlines1)
        resHTMLlines1 = re.sub('<td colspan="2"><b>1st tree:</b> Orthogroup</td>', '<td colspan="2"><b>1st tree with mode E1st:</b> Speciation/duplication events in the query sequnece lineage</td>', resHTMLlines1)
    if mode in ("ReE1st", "D1st") and BSthreshold_4_ReE1st == "reconcile":
        print("13121212")
        resHTMLlines1 = re.sub('<td align="center" valign="top">Rearranged gene tree.*</td>', '<td align="center" valign="top">NONE</td>', resHTMLlines1)
        resHTMLlines1 = re.sub('.*2ndGeneTree.pdf.*', '<td align="center" valign="top">NONE</td>', resHTMLlines1)
        resHTMLlines1 = re.sub('<td align="center" valign="top" name="REARRANGED2">Rearranged gene tree', '<td align="center" valign="top" name="REARRANGED2">Reconciled gene tree', resHTMLlines1)

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

    #out = open(eachDirAddress_FN + "300_resultsREA.html", "w")
    out = open(queryID + ".html", "w")
    out.write(resHTMLlines1)
    out.close()

### result_HTML End
###############################################################


###############################################################
### Others START

def deleteFiles(eachDirAddress):

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
            if re.search(r"^" + keyword, fileName):
                fileNames_rm.append(fileName)
    for file in fileNames_rm:
        address_file = eachDirAddress + file
        line_rm = "rm " + address_file
        #print("line_rm", line_rm)
        subprocess.call(line_rm, shell=True)

### Others END


########### Data summarize after all gene tree estimated
def make_lines_atmarkSeparated(
    mode,
    SpeciesTree,
    keyNode,
    allNodes_speciesTree,
    childSpeciesNodes_AllGroup,
    childSpeciesNodes_focalGroup,
    speciesWithGeneFunction,
    outdir, geneIDs_fn
):
    print("### make_lines_atmarkSeparated() ###")

    speciesNames_in_orthogroup = collect_speciesNames_in_orthogroup(SpeciesTree, keyNode, allNodes_speciesTree)


    #querySpeciesNode = identify_speciesNode(name_querySpecies)
    #speciesNodes_including_querySpecies = collect_ancestralNodes(allNodes_speciesTree, querySpeciesNode)

    lines_FN = []
    for x in range(len(geneIDs_fn)):
        print(x+1, geneIDs_fn[x])
        
        line_FN = "Num" + "@" + str(x+1) + ","

        line_FN += "QueryGeneID" + "@" +geneIDs_fn[x] + ","

        address_eachDirectory = os.path.join(outdir, geneIDs_fn[x])
        if not os.path.exists(address_eachDirectory):
            print("Cannot find directories. Stopped:")
            print(address_eachDirectory)
            exit()

        address_100_resultsREA_html = os.path.join(outdir, geneIDs_fn[x], "100_analysisSummary.txt")
        if not os.path.exists(address_100_resultsREA_html):
            print("Cannot find 100_analysisSummary.txt. Stopped:")
            print(address_100_resultsREA_html)
            exit()

        address_100_analysisSummary_txt = os.path.join(outdir, geneIDs_fn[x], "100_analysisSummary.txt")
        if not os.path.exists(address_100_analysisSummary_txt):
            print("Cannot find 100_analysisSummary.txt. Stopped:")
            print(address_100_analysisSummary_txt)
            exit()

        seqDict = readRes_dict(address_100_analysisSummary_txt)

        #for name, val in seqDict.items():
        #    print("name:", name)
        #    #print("val:", val)
        #exit()
    
        line_FN += "BS_of_orthogroupBasalNode@"
        if ">BS_of_orthogroupBasalNode" in seqDict.keys():
            line_FN += seqDict[">BS_of_orthogroupBasalNode"][0] + ","
        else:
            line_FN += "NONE" + ","

        line_FN += "QueryLength@"
        if ">NumberAssigned_querySequence" in  seqDict.keys():
            line_FN += str(len(seqDict[">NumberAssigned_querySequence"][1])) + ","
        else:
            line_FN += "NONE" + ","

        line_FN += "SpeciesWithGeneFunction@"
        if ">Orthogroup" in seqDict.keys():
            flagTMP = 0
            for geneLeaf in seqDict[">Orthogroup"]:
                if re.search(r"^" + speciesWithGeneFunction + "_", geneLeaf):
                    flagTMP += 1
                    line_FN += geneLeaf + " "
            if flagTMP < 1:
                line_FN += "NONE,"
            else:
                line_FN += ","
            #exit()
        else:
            line_FN += "NONE" + ","
    
        line_FN += "2ndGeneTree@"
        if ">2nd_rearranged_gene_tree_newick" in seqDict.keys():
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

        
        if ">Number_of_blastHits" in seqDict.keys():
            for node_num in seqDict[">Number_of_blastHits"]:
                match = re.search(r"^([^ ]+) +(\d+)$", node_num)
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

        if ">GeneNumber_of_orthogroup" in seqDict.keys():
            for node_num in seqDict[">GeneNumber_of_orthogroup"]:
                match = re.search(r"^([^ ]+) +(\d+)$", node_num)
                node = match.group(1)
                num = match.group(2)
                line_FN += "OGnum_" + node + "@" + num + ","
        else:
            for speciesName_in_orthogroup in speciesNames_in_orthogroup:
                #line_FN += "OGnum_" + speciesName_in_orthogroup + "@NONE" + ","
                #line_FN += "OGnum_" + speciesName_in_orthogroup + "@" + " " + ","
                line_FN += "OGnum_" + speciesName_in_orthogroup + "@NONE" + ","


        if mode == "S":
            key_nameLine_monophyly = ">MonophyleticGeneGroups"
            key_nameLine_sister = ">SisterGeneGroups"
            childSpeciesNodes = childSpeciesNodes_focalGroup
        else:
            key_nameLine_monophyly = ">MonophyleticGeneGroups_1stTree"
            key_nameLine_sister = ">SisterGeneGroups_1stTree"
            childSpeciesNodes = childSpeciesNodes_AllGroup

        if key_nameLine_monophyly in seqDict.keys():
            #print("seqDict[>Number_of_duplicatedNode]", seqDict[">Number_of_duplicatedNode"])
            for node_BS_sister in seqDict[key_nameLine_monophyly]:
                #print("node_BS_sister", node_BS_sister, "|")
                match = re.search(r"^([^ ]+) +([^ ]+) +([^ ]+)$", node_BS_sister)   ##################
                #match = re.search(r"^([^ :]+)[ :]+([^ ]+)$", node_BS_sister)       ##################
                name_targetNode = match.group(1)
                bsBaclue = match.group(2)
                duplicationStatus = match.group(3)
                line_FN += "BS_of_" + name_targetNode + "_monophyly@" + bsBaclue + ","
                line_FN += "dupStatus_" + name_targetNode + "@" + duplicationStatus + ","
        else:
            #for targetSpeciesNode in speciesNodes_including_querySpecies:
            #for targetSpeciesNode in childSpeciesNodes_focalGroup:
            for targetSpeciesNode in childSpeciesNodes:
                name_targetSpeciesNode = make_nodeName_from_nodeLavel_NHXstyle(targetSpeciesNode[2])
                #print("name_targetSpeciesNode", name_targetSpeciesNode)
                line_FN += "BS_of_" + name_targetSpeciesNode + "_monophyly@NONE" + ","
                line_FN += "dupStatus_" + name_targetSpeciesNode + "@NONE" + ","

        if key_nameLine_sister in seqDict.keys():
            #print("seqDict[>Number_of_duplicatedNode]", seqDict[">Number_of_duplicatedNode"])
            for node_BS_sister in seqDict[key_nameLine_sister]:
                #print("node_BS_sister", node_BS_sister, "|")
                match = re.search(r"^([^ ]+) +([^ ]+) +([^ ]+)$", node_BS_sister)   ##################
                #match = re.search(r"^([^ :]+)[ :]+([^ ]+) +([^ ]+)$", node_BS_sister)   ##################
                name_targetNode = match.group(1)
                bsBaclue = match.group(2)
                name_sisterNode = match.group(3)
                line_FN += "Sister_of_" + name_targetNode + "@" + name_sisterNode + ","
                line_FN += "BS_with_" + name_targetNode + "@" + bsBaclue + ","
        else:
            #for targetSpeciesNode in speciesNodes_including_querySpecies:
            for targetSpeciesNode in childSpeciesNodes:
                name_targetSpeciesNode = make_nodeName_from_nodeLavel_NHXstyle(targetSpeciesNode[2])
                #print("name_targetSpeciesNode", name_targetSpeciesNode)
                line_FN += "Sister_of_" + name_targetSpeciesNode + "@NONE" + ","
                line_FN += "BS_with_" + name_targetSpeciesNode + "@NONE" + ","

        #if ">Number_of_duplicatedNode" in  seqDict.keys():
        #    #print("seqDict[>Number_of_duplicatedNode]", seqDict[">Number_of_duplicatedNode"])
        #    for node_dup in seqDict[">Number_of_duplicatedNode"]:
        #        #print("node_dup", node_dup)
        #        match = re.search(r"^([^ ]+) +(\d+)$", node_dup)
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
            line = re.sub(r",[^@]+@", ",", line)   ########################### comment out for make sure
            line = re.sub(r"^[^@]+@", "", line)    ########################### comment out for make sure
            #print("line2", line)
        out.write(line + "\n")
    out.close()


def make_indexLine(first_line):
    indexLine_FN = re.sub(r"@[^,]+,", ",", first_line)
    indexLine_FN = re.sub(r"@[^@]+$", "", indexLine_FN)
    indexLine_FN = re.sub(r"@,", ",", indexLine_FN)
    indexLine_FN = re.sub(r"@$", "", indexLine_FN)
    return indexLine_FN

def copy_alignment_orthogroup(eachDirAddress, dataset, queryID, alignment_orthogroups):
    if dataset == "Exclude3rd" or dataset == "Include3rd":
        if os.path.exists(eachDirAddress + "180_aln_nucl_fas.txt"):
            #print("Present")
            recsTMP = readFasta_dict(eachDirAddress, "180_aln_nucl_fas.txt")
            #print("queryID", queryID)
            out = open(alignment_orthogroups + "/" + queryID + ".txt", "w")
            for nameLine, seq in recsTMP.items():
                out.write(nameLine + "\n")
                out.write(seq + "\n")
            out.close()
    else:
        if os.path.exists(eachDirAddress + "190_aln_prot_fas.txt"):
            #print("Present")
            recsTMP = readFasta_dict(eachDirAddress, "190_aln_prot_fas.txt")
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
def main():
    args = parse_args(sys.argv[1:])
    query_id = args.query_id

    ctx = initialize_context(query_id)
    #print("ctx.eachDirAddress_e1stre", ctx.eachDirAddress_e1stre)
    #print("### line 4550")
    #exit()

    echo_major_settings(ctx)
    prepare_files_and_species_tree(ctx)

    if ctx.mode in ("S", "S1st", "ReS1st"):
        run_mode_S_block(ctx)     # ここで終了

    if ctx.mode in ("D", "D1st", "ReD1st"):
        run_mode_D_draw_only(ctx) # ここで終了

    #print("ctx.eachDirAddress", ctx.eachDirAddress)
    #print("### line 4353")
    #exit()
    if ctx.mode in ("E", "E1st", "ReE1st"):
        run_mode_E_pipeline(ctx)      # E/E1st/ReE1st

if __name__ == "__main__":
    main()

exit()
