import numpy
import matplotlib.pyplot as plt
import subprocess

meme_matrix_template = "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\nMOTIF\t%s\nletter-probability matrix: alength= 4 w= 18 nsites= 1 E= 0\n%s"

def run_cmd(cmd):
    subprocess.call(cmd, shell=True)


def select_opt_motif(report_path):

    with open(report_path) as report:
        header = next(report) # skip and store header
        # 0:EXPERIMENT,1:MOTIF,2:OPTIMIZED,3:INITIAL_AUC,
        # 4:INITIAL_P,5:OPTIMIZED_AUC,6:OPTIMIZED_P,7:CORRELATION,8:CORRELATION_P
        corr_p = int(8)
        opt_p = int(6)
        # Mock motif
        selected_motif = ['', 'error', '0', '0.5', '1', \
                          '-1', '1', '0', '1']

        for line in report:
            current_motif = line.replace("\n", "").split("\t")

            # Do not inclde motifs without correlation to their intial sequence.
            if float(current_motif[corr_p]) > float(0.001):
                continue

            # Redefine selected motif as the one with the smallest pval.
            if float(current_motif[opt_p]) < float(selected_motif[opt_p]):
                selected_motif = current_motif


        if selected_motif[1] == "error":
            print("There is no motif with a CORRELATION_P smaller than 0.001.\nStopping")
            exit()
            
    return selected_motif

def convert_cisbp_to_meme(matrix):

    name = "None"
    position = False
    frequencies = ""
    matrix_ = matrix.split("\n")

    for line in matrix_:
        
        if line.find("Motif\t") >= 0:
            name = line.split("\t")[1]
            continue
            
        if position:
            line_ = "\t".join(line.split("\t")[1:])            
            frequencies += "%s\n" % (line_)

        if line == "Pos\tA\tC\tG\tT":
            position = True

    return meme_matrix_template % (name, frequencies)


def write_cisbp_meme_motif(motif_id, motif_cisbp_file, motif_meme_file, PFMS_PATH):

    with open(PFMS_PATH, 'r') as matrices:
        matrices = matrices.read().split("\n\n")
        motif_id_line = "Motif\t" + motif_id + "\n"
        
        for matrix in matrices:
            
            if matrix.find(motif_id_line) >= 0:
                
                cisbp_motif = open(motif_cisbp_file, 'w')
                cisbp_motif.write(matrix)
                cisbp_motif.close()

                meme_motif = open(motif_meme_file, 'w')
                meme_matrix = convert_cisbp_to_meme(matrix)
                meme_motif.write(meme_matrix)
                meme_motif.close()



def smooth(y, box_pts):
    box = numpy.ones(box_pts)/box_pts
    y_smooth = numpy.convolve(y, box, mode='same')
    return y_smooth

def write_centrimo_plot(site_table, out_centrimo_graph, centrimo_score):
    points, sites, counts, norm_counts, sum_counts = [], [], [], [], 0.0
    
    with open(site_table, 'r') as centrimo_table:
        next(centrimo_table) # Header
        for line in centrimo_table:
            line = line.replace("\n", "").split("\t")
            sites.append(float(line[0]))
            counts.append(float(line[1]))
            sum_counts += float(line[1])
            
    for count in counts:
        # Normalization
        norm_counts.append(count / sum_counts)
        
    fig, ax = plt.subplots(figsize=(10,5), facecolor=("white"))
    ax.set_facecolor("lightgray")
    ax.set_title("Motif Probability Graph (score ≥ %s bits)" % centrimo_score)
    # ax.plot(sites, smooth(norm_counts, 3), color='red', lw=0.2)
    ax.plot(sites, smooth(norm_counts, 20), color='navy', lw=1.5)
    ax.set_xlabel('Position of Best Site In Sequence')
    ax.set_ylabel('Probability')
    plt.grid()
    
    ## Supported formats: png, pdf, ps, eps and svg
    plt.savefig(out_centrimo_graph, dpi=100, facecolor='White', edgecolor='w', \
                orientation='portrait', format="png", pad_inches=0.01) 
        

def write_html(metadata, centrimo_graph, PREFIX, job_id):
    # ['', 'CTCF_HUMAN:3-8', '1', '0.814348', '1.20937e-66',
    #  '0.867272', '3.51351e-90', '0.86303', '1.55052e-07']
    # 0:EXPERIMENT,1:MOTIF,2:OPTIMIZED,3:INITIAL_AUC,
    # 4:INITIAL_P,5:OPTIMIZED_AUC,6:OPTIMIZED_P,7:CORRELATION,8:CORRELATION_P
    motif_name = metadata[1].replace("-", "_").replace(":", "_")

    ## Logo PNG.
    logo_png = "logo" + motif_name + ".png"
    logo_opt_png = "logo" + motif_name + "_opt.png"
    ## SEED
    seed_auc = metadata[3]
    seed_pval = metadata[4]
    ## Optimized
    opt_auc = metadata[5]
    opt_pval = metadata[6]
    ## CORR
    pcorr = metadata[7]
    corr_pval = metadata[8]
    
    # job_id, motif_name, seed_auc, seed_pval, opt_auc, opt_pval
    # pearson_cor, cor_pval, seed_png, opt_png, sites_png
    template = load_html_template()
    
    html_file = template % (job_id, job_id, motif_name, seed_auc, \
                            seed_pval, opt_auc, \
                            opt_pval, pcorr, corr_pval, logo_png, \
                            logo_opt_png, centrimo_graph)

    html_path = PREFIX + "/report.html"
    OUT_html = open(html_path, 'w')
    OUT_html.write(html_file)
    OUT_html.close()

    
def load_html_template():

    template = """
    <HEAD>
    <TITLE>RCADE2 %s</TITLE>
    </HEAD>
    <BODY BGCOLOR='WHITE'>
    <CENTER>
    <H1>RCADE2 <BR> %s</H1>
    </CENTER>
    <head>
    <style>
    table {
    font-family: arial, sans-serif;
    border-collapse: collapse;
    width: 100%%;
    }
    
    td, th {
    border: 1px solid #dddddd;
    text-align: left;
    padding: 8px;
    }
    tr:nth-child(even) {
    background-color: #dddddd;
    }
    </style>
    </head>
    <body>
    <H2>%s</H2>
    <table>
    <tr>
    <th></th>
    <th>AUC</th>
    <th>p-value</th>
    </tr>
    <tr>
    <td>Seed</td>
    <td>%s</td>
    <td>%s</td>
    </tr>
    <tr>
    <td>Optimized</td>
    <td>%s</td>
    <td>%s</td>
    </tr>
    </table>
    <H4>Pearson correlation: %s, p-value: %s</H4>
    <CENTER>
    <H2>Seed</H2>
    <IMG SRC='%s'>
    <H2>Optimized</H2>
    <IMG SRC='%s'>
    </CENTER>
    <BR>
    <H1>Site counts</H1>
    <CENTER>
    <IMG SRC='%s'>
    </CENTER>
    <BR>
    <H4>References:</H4>
    Najafabadi, H. S., Albu, M., & Hughes, T. R. (2015). Identification of C2H2-ZF binding preferences from ChIP-seq data using RCADE. Bioinformatics, 31(17), 2879-2881. doi:10.1093/bioinformatics/btv284. <a href="https://academic.oup.com/bioinformatics/article/31/17/2879/183644">[full text]</a><BR>  
    Bailey, T. L., & Machanick, P. (2012). Inferring direct DNA binding from ChIP-seq. Nucleic Acids Res, 40(17), e128. doi:10.1093/nar/gks433 <a href="http://nar.oxfordjournals.org/content/40/17/e128">[full text]</a>        
    </body>
    """
    return template

