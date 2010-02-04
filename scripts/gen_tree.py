#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import optparse
import shutil
from ete2 import Tree, faces
sys.path.append('/users/rg/agrimaldi/usr/lib/python2.5/site-packages')
from AGBio.selenoprofiles_tools.results_analyser import *

MAX_SN_LEN = 10
MAX_LN_LEN = 30

species = {'Candidatus_Solibacter_usitatus_Ellin6076_' : 'Solibacter_usitatus_Ellin6076_',
           'Candidatus_Koribacter_versatilis_Ellin345_' : 'Acidobacteria_bacterium_Ellin345_',
           'Synechococcus_sp._JA-2-3Ba_2-13_' : 'Synechococcus_sp._JA-2-3Ba2-13_',
           'Baumannia_cicadellinicola_str._Hc__Homalodisca_coagulata_' : 'Baumannia_cicadellinicola_str._Hc_Homalodisca_coagulata_',
           'Cronobacter_sakazakii_ATCC_BAA-894_' : 'Enterobacter_sakazakii_ATCC_BAA-894_',
           'Mycoplasma_gallisepticum_str._R_' : 'Mycoplasma_gallisepticum_R_',
           'Leptospira_biflexa_serovar_Patoc_strain_Patoc_1__Paris_' : 'Leptospira_biflexa_serovar_Patoc_strain_Patoc_1_Paris_',
           'Leptospira_biflexa_serovar_Patoc_strain_Patoc_1__Ames_' : 'Leptospira_biflexa_serovar_Patoc_strain_Patoc_1_Ames_',
           'Chlorobium_luteolum_DSM_273_' : 'Pelodictyon_luteolum_DSM_273_',
           'Chelativorans_sp._BNC1_' : 'Mesorhizobium_sp._BNC1_',
           'Candidatus_Ruthia_magnifica_str._Cm__Calyptogena_magnifica_' : 'Candidatus_Ruthia_magnifica_str._Cm_Calyptogena_magnifica_',
           'Alkalilimnicola_ehrlichii_MLHE-1_' : 'Alkalilimnicola_ehrlichei_MLHE-1_',
           'Escherichia_coli_O157_H7_EDL933_' : 'Escherichia_coli_O157-H7_EDL933_',
           'Escherichia_coli_O157_H7_str._Sakai_' : 'Escherichia_coli_O157-H7_str._Sakai_',
           'Salmonella_enterica_subsp._arizonae_serovar_62_z4_z23_' : 'Salmonella_enterica_subsp._arizonae_serovar_62-z4-z23_',
           'Buchnera_aphidicola_str._Bp__Baizongia_pistaciae_' : 'Buchnera_aphidicola_str._Bp_Baizongia_pistaciae_',
           'Buchnera_aphidicola_str._APS__Acyrthosiphon_pisum_' : 'Buchnera_aphidicola_str._APS_Acyrthosiphon_pisum_',
           'Buchnera_aphidicola_str._Sg__Schizaphis_graminum_' : 'Buchnera_aphidicola_str._Sg_Schizaphis_graminum_',
           'Buchnera_aphidicola_str._Cc__Cinara_cedri_' : 'Buchnera_aphidicola_str._Cc_Cinara_cedri_',
           'Streptomyces_coelicolor_A3_2_' : 'Streptomyces_coelicolor_A32_',
           'Polynucleobacter_necessarius_subsp._asymbioticus_QLW-P1DMWA-1_' : 'Polynucleobacter_sp._QLW-P1DMWA-1_',
           'Polynucleobacter_necessarius_subsp._necessarius_STIR1_' : 'Polynucleobacter_necessarius_STIR1_',
           'Acidovorax_citrulli_AAC00-1_' : 'Acidovorax_avenae_subsp._citrulli_AAC00-1_',
           'Ruegeria_sp._TM1040_' : 'Silicibacter_sp._TM1040_',
           'Ruegeria_pomeroyi_DSS-3_' : 'Silicibacter_pomeroyi_DSS-3_',
           'Brucella_melitensis_bv._1_str._16M_' : 'Brucella_melitensis_16M_',
           'Yersinia_pestis_KIM_10_' : 'Yersinia_pestis_KIM_',
           'Escherichia_coli_str._K-12_substr._DH10B_' : 'Escherichia_coli_str._K12_substr._DH10B_',
           'Escherichia_coli_str._K-12_substr._W3110_' : 'Escherichia_coli_str._K12_substr._W3110_',
           'Escherichia_coli_str._K-12_substr._MG1655_' : 'Escherichia_coli_str._K12_substr._MG1655_',
           'Salmonella_enterica_subsp._enterica_serovar_Typhimurium_str._LT2_' : '',
           'Rhodococcus_jostii_RHA1_' : 'Rhodococcus_sp._RHA1_',
           'Bacillus_cytotoxicus_NVH_391-98_' : 'Bacillus_cereus_subsp._cytotoxis_NVH_391-98_',
           'Pseudomonas_fluorescens_Pf0-1_' : 'Pseudomonas_fluorescens_PfO-1_',
           'Desulfovibrio_vulgaris_str._Hildenborough_' : 'Desulfovibrio_vulgaris_subsp._vulgaris_str._Hildenborough_',
           'Desulfovibrio_vulgaris_DP4_' : 'Desulfovibrio_vulgaris_subsp._vulgaris_DP4_'
           }

def long2short(longname):
    ln = longname.split()
    return ln[0][0] + ln[1][:3]

def add_to_species_dict(longname, genomes):
    sln = sanitize(longname)
    if sln in genomes and sln not in species:
        species[sln] = sln

def sanitize(longname):
    t = longname
    for c in '\\\'\"/':
        t = t.replace(c, '')
    for c in '-_':
        t = t.rstrip(c)
    t = t.replace(' ', '_')
    return t + '_'

def layout(node):
    if node.is_leaf():

        add_to_species_dict(node.name, g_genomes)
        node.img_style['size'] = 10
        shortNameFace = faces.TextFace(long2short(node.name).ljust(MAX_SN_LEN), ftype='monospace')
        pathNameFace = faces.TextFace(species[sanitize(node.name)], ftype='monospace')
        faces.add_face_to_node(shortNameFace, node, column=0, aligned=True)
        faces.add_face_to_node(pathNameFace, node, column=1, aligned=True)

        try:
            fp = os.path.join(resfolder, species[sanitize(node.name)])
            sp_parser = GenomeFolderParser(fp)
            sp_parser.parse(sec=True, cys=True)
            if sp_parser.cys and not sp_parser.sec:
                node.img_style['bgcolor'] = '#9db0cf'
            if sp_parser.sec:
                node.img_style['bgcolor'] = '#ff503c'
        except Exception:
            pass

    else:
        node.img_style['size'] = 0
        

def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--input_file',
                      dest='infilename',
                      help='input filename, newick tree.',
                      metavar='FILE')
    
    parser.add_option('-s', '--selenoprofiles_results_folder',
                      dest='sp_res_folder',
                      help='result folder of a selenoprofiles run.',
                      metavar='DIR')

    parser.add_option('-g', '--gui',
                      action='store_true', dest='gui', default=False,
                      help='Use GUI for printing.')
    
    (options, args) = parser.parse_args()

    global g_genomes
    global resfolder

    g_genomes = os.listdir(options.sp_res_folder)

    resfolder = options.sp_res_folder
    
    t = Tree(options.infilename)

    if options.gui:
        t.show(layout)
    else:
        print t


if __name__ == '__main__':
    try:
        main()
    except Exception, (e):
        print e
        
