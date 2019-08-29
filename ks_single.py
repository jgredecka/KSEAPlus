def userInput(ks_db, graphics, df, min_sub):

    import io
    from io import StringIO
    import pandas as pd
    import scipy.stats as st
    import numpy as np
    
    if graphics == "yes":
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt

    # Function that returns the Kolmogorov-Smirnov test statistic along with the p-value.
    # -log10 of p-value is also calculated for the barplot.
    def ksTest(substrates, non_substrates):
        result = st.ks_2samp(substrates, non_substrates)
        ks_stat = result[0]
        pval = result[1]
        log_pval = np.log10(1/pval)
        return ks_stat, pval, log_pval

    # Convert JSON-string datatset back into a DF.
    df = pd.read_json(df, orient="split")
    
    # User data is passed from the server and parsed as appropriate.
    user_file=df.values.tolist()
    array=[]
    for line in user_file:
        phosphosite=line[0]
        fc_value = line[1]
        if phosphosite == "" or fc_value == "":
            continue
        else:
            array.append([phosphosite, float(fc_value)])

    # Multiple phosphosites separated by a colon are split here.
    # This ensures each phosphosite substrate starts with a new line.
    data=[]
    for n in range (0, len(array)):
        site=array[n][0].upper()
        fc=array[n][1]
        site=site.split(";")
        for s in site:
            if s == '':
                continue
            else:
                data.append([s, fc])

    # Mapping of phosphosites to their (often multiple) log2(FC) values is achieved here.
    dic={}
    for entry in data:
        site=entry[0]
        fc=entry[1]
        if site not in dic:
            dic[site]=[fc]
        else:
            dic[site].append(fc)

    # If the same phosphosite has been detected more than once, a mean log2(FC) is calculated for that phosphosite.
    # Final array contains unique phosphosites and individual log2(FC) values, averaged where appropriate.
    for key in dic:
        length=len(dic[key])
        mean_fc=sum(dic[key])/length
        dic[key] = float(mean_fc)

    # Each phosphosite in dic is scanned against the K-S db. 
    # If a match is found, relevant information for that phosphosite is appended to ks_links.
    ks_links=[]
    for x in dic:
        for y in ks_db:
            if x == y[0]:
                ks_links.append([y[1], y[0], y[2], y[3], dic[x]])
            else:
                continue

    # The array is then converted into a dataframe to be viewed as a table.
    ks_links_df = pd.DataFrame(ks_links, columns=["Kinase", "Site", "Site.Seq(+/- 7AA)", "Source", "log2(FC)"])

    # A dictionary containing unique kinases as keys and substrate log2(FCs) as values is created.
    # If the same kinase was identified for multiple substrates, multiple FCs are appended to the dictionary values.
    kinase_dic={}
    for match in ks_links:
        kinase=match[0]
        log2fc=match[4]
        if kinase not in kinase_dic:
            kinase_dic[kinase]=[log2fc]
        else:
            kinase_dic[kinase].append(log2fc)

    # The dictionary is used to calculate the number of substrates identified for each unique kinase.
    # It also calculates the mean log2(FC) across each kinase's substrates.
    # The algorithm computes a signed KS statistic, p-value and -log10(p-value) using substrate and non-substrate log2(FC) values.    
    kol_smir_info=[]
    for kinase in kinase_dic:
        substrate_num=len(kinase_dic[kinase])
        kin_fc_mean=sum(kinase_dic[kinase]) / float(substrate_num)
        # Substrate names corresponding to a given kinase are extracted into a flat list.
        # These will be used to identify all non-substrates conditionally. 
        sub_df = ks_links_df.loc[ks_links_df['Kinase'] == kinase]
        sub_names = sub_df['Site'].tolist()
        # log2(FC) values for substrates.
        sub_fc = kinase_dic[kinase]
        # Non-substrates are located within the ks_links dataframe.
        # A sub-dataframe containing non-substrates only is generated.
        # log2(FC) values associated with the non-substrates are then stored in a new array.
        non_sub_df=ks_links_df.loc[~ks_links_df['Site'].isin(sub_names)]
        non_sub_df=non_sub_df.drop_duplicates(subset='Site')
        non_sub_fc = non_sub_df['log2(FC)'].tolist()
        # KS-test statistic and p-value for each kinase are calculated here for kinases with at least 3 substrates.
        ks_stat, pval, log_pval = ksTest(sub_fc, non_sub_fc)
        # -log10(p-val) and KS is signed based on the mean log2(FC) of the kinase.
        if kin_fc_mean < 0:
            log_pval = -log_pval
            ks_stat = -ks_stat
        # All relevant results are stored in a new array and converted into a dataframe.
        kol_smir_info.append([kinase, substrate_num, kin_fc_mean, ks_stat, pval, log_pval])

    # KSEA results are converted into a dataframe.
    kol_smir_df = pd.DataFrame(kol_smir_info, columns = ['Kinase', 'Sub.Count', 'mnlog2(FC)', '(+/-) KS', 'pVal', "(+/-) -log10(pVal)"])

    # Barplot is only generated if the user chose to produce graphics at file upload.
    if graphics == "no":
        svg_fig = "Barplot was not generated for this analysis."
    elif graphics == "yes":
        # The array is sorted according to the -log10(p-value) to improve plot readability.
        kol_smir_info.sort(key=lambda x: float(x[5]))

        # Extraction of desirable results that match user-defined parameters.
        # min_sub denotes the minimum substrate count per kinase indicated by the user and must be >= 3.
        kinases=[]
        log_pvals=[]
        for entry in kol_smir_info:
            kin=entry[0]
            sub=entry[1]
            logp = entry[5]
            if sub >= min_sub:
                kinases.append(kin)
                log_pvals.append(logp)

        # Dictionary of kinases as keys and p-values as values. 
        # It will serve as a map for the colour coding of the barplot bars.
        kin_map={}
        for entry in kol_smir_info:
            kin=entry[0]
            pval=entry[4]
            kin_map[kin]=pval

        # A list of colours for the barplot. 
        colours=[]
        for kin in kinases:
            # Bars are pale red if p-value is smaller than 0.05.
            if kin_map[kin] < 0.05 and kin_map[kin] >= 0.01:
                colours.append("#ffb3b3")
            # Bars are bright red if p-value is less than 0.01.
            elif kin_map[kin] < 0.01:
                colours.append("#ff4d4d")
            # Otherwise the bars remain grey.
            else:
                 colours.append("#b3b3b3")

        # Generate an array of sequence to be used to plot the y-axis values.
        y_pos = np.arange(len(kinases))
        n=len(y_pos)

        # Set the margins and bar height for a single category.
        topmargin = 0.1 #inches
        bottommargin = 0.1 #inches
        categorysize = 0.15 # inches

        # Calculate a dynamic figure height based on the known values above.
        figheight = topmargin + bottommargin + (n+1)*categorysize

        # Figure box and axes are plotted here.
        fig, ax = plt.subplots(figsize=(5, figheight))

        # Horizontal bars are generated.
        ax.barh(y_pos, log_pvals, color=colours)

        # Y-axis limit to reduce the top and bottom margin white spaces.
        plt.ylim(min(y_pos)-1, max(y_pos)+1)

        # Barplot is styled and labelled.
        ax.set_yticks(y_pos)
        ax.set_yticklabels(kinases)
        for pos in ['right','top', 'left']:
            ax.spines[pos].set_visible(False)
        ax.tick_params(left=False)
        plt.yticks(fontsize=6)
        plt.xticks(fontsize=6)
        plt.ylabel("Kinase", fontsize=6)
        plt.xlabel("(+/-) -log10(p-value)", fontsize=6)

        # A legend to indicate the meaning of the bar colours.
        non_sig = mpatches.Patch(color="#b3b3b3", label='p-value >= 0.05')
        sig = mpatches.Patch(color="#ffb3b3", label='p-value >= 0.01')
        extra_sig = mpatches.Patch(color="#ff4d4d", label='p-value < 0.01')
        plt.legend(handles=[non_sig, sig, extra_sig], bbox_to_anchor=(-0.45, 1), loc=2, fontsize=6)

        # Create a StringIO object and use it to write SVG figure data to string buffer.
        fig_file = StringIO()
        fig.savefig(fig_file, format='svg', bbox_inches="tight")
        # Seek beginning of the figure file.
        fig_file.seek(0)
        # Retrieve figure contents as a string.
        svg_fig = '<svg' + fig_file.getvalue().split('<svg')[1]
        # Free memory buffer.
        fig_file.close()

    # Convert results DFs into JSON strings.
    kol_smir_df = kol_smir_df.to_json(orient='split')
    ks_links_df = ks_links_df.to_json(orient='split')
    
    return kol_smir_df, ks_links_df, svg_fig