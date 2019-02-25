#!/bin/env python3
import os
import argparse
import my_utils

parser = argparse.ArgumentParser(description='Create a HTML summary page for RCADE2')

parser.add_argument('--peaks', dest='PEAK_PATH', action="store", \
                    type=str, \
                    help='Path to peak file.')

parser.add_argument('--report', dest='REPORT_PATH', action="store",\
                    type=str, \
                    help='PFMs report table of tested motifs.')

parser.add_argument('--pfms', dest='PFMS_PATH', action="store", \
                    type=str,
                    help='Path to PFMs, all seeds/optimized motifs in CIS-BP format.')

parser.add_argument('--prefix', dest='PREFIX', action="store", \
                    type=str,
                    help='Out path prefix for results.')

parser.add_argument('--job_id', dest='JOB_ID', action="store", \
                    type=str, default="JOB_ID", \
                    help='Job ID for RCADE2.')

args = parser.parse_args()

## Header example for results..report.txt
# HEADER=["EXPERIMENT", "MOTIF", "OPTIMIZED", "INITIAL_AUC", \
#         "INITIAL_P", "OPTIMIZED_AUC", "OPTIMIZED_P", \
#         "CORRELATION", "CORRELATION_P"]

def main():

    if os.path.getsize(args.REPORT_PATH) == 0 or \
       os.path.getsize(args.PFMS_PATH) == 0 or \
       os.path.getsize(args.PEAK_PATH) == 0:
        print("Input files are empty.")
        exit()
       
        
    log_file = open(args.PREFIX + "/centrimo/centrimo.log", "w")
    ###########################################################################
    ## Return the information of the selected motif, name and scores
    selected_motif = my_utils.select_opt_motif(args.REPORT_PATH)

    ## Base name of the selected motif.
    seed_id = str(selected_motif[1])
    ## Name of the optimized selected motif.
    opt_id = seed_id + "|opt"

    ###########################################################################
    ## Selected seed  motif                                                  ##
    ###########################################################################
    # Out file name: results.seed.PFM.txt
    seed_cisbp_path = args.PREFIX + "/results.seed.PFM.txt"
    # Out file name: results.opt.PFM.meme.txt
    seed_meme_path = args.PREFIX + "/results.seed.PFM.meme.txt"

    ## Use name to write selected motif to one file.
    my_utils.write_cisbp_meme_motif(seed_id, seed_cisbp_path, seed_meme_path, \
                                    args.PFMS_PATH)
    
    ###########################################################################
    ## Selected optimized motif                                              ##
    ###########################################################################
    # Out file name: results.seed.PFM.txt
    opt_cisbp_path = args.PREFIX + "/results.opt.PFM.txt"
    ## Convert matrix to meme format
    opt_meme_path = args.PREFIX + "/results.opt.PFM.meme.txt"
    
    ## Use name to write selected motif in two formats/files
    my_utils.write_cisbp_meme_motif(opt_id, opt_cisbp_path, opt_meme_path, \
                                    args.PFMS_PATH)
    
    ###########################################################################
    ## Run Centrimo                                                          ##
    ###########################################################################
    centrimo_score = 5 # Default.
    # Must be directory as it creates another directory in it.
    centrimo_out = args.PREFIX + "/centrimo"
    cmd_centrimo = "centrimo --oc %s --score %s %s %s &>> log.step3.txt" % \
                   (centrimo_out, centrimo_score, args.PEAK_PATH, opt_meme_path)

    print("Running centrimo", file=log_file)
    my_utils.run_cmd(cmd_centrimo)
    site_table = centrimo_out + "/site_counts.txt"
    centrimo_graph_name = "motif_sites_graph.png"
    out_centrimo_graph = args.PREFIX + "/" + centrimo_graph_name
    my_utils.write_centrimo_plot(site_table, out_centrimo_graph, centrimo_score)
    
    ###########################################################################
    ## Create PNG image of the logo for optimized and normal motifs          ##
    ###########################################################################
    ####  Selected seed
    cmd_seed_image = "meme2images -png %s %s" % (seed_meme_path, args.PREFIX)
    my_utils.run_cmd(cmd_seed_image)

    ####  Selected optimized
    cmd_opt_image = "meme2images -png %s %s" % (opt_meme_path, args.PREFIX)
    my_utils.run_cmd(cmd_opt_image)

    print("Creating HTML report", file=log_file)
    ############################################################
    my_utils.write_html(selected_motif, centrimo_graph_name, args.PREFIX, args.JOB_ID)
    ############################################################
    rm_cmd = "rm %s/results.opt.ps" % (args.PREFIX)
    my_utils.run_cmd(rm_cmd)
    
if __name__ == "__main__":

    main()
