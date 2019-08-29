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

    # Function used to convert kinase z-scores to corresponding p-values.
    def getpValue(z):
        if z < 0:
            p=st.norm.cdf(z)
            return p
        else:
            dist=st.norm.cdf(z)
            p=1.0 - dist
            return p
      
    # Function to calculate the Z-score for each kinase.
    def zScore(mean_kin, mean_all, sub_num, sd):
        z = (mean_kin - mean_all) * sub_num**(1/2) / sd
        return z
    
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

    # Mapping of phosphosites to their (often multiple) log2(FCs) is achieved here.
    dic={}
    for entry in data:
        site=entry[0]
        fc=entry[1]
        if site not in dic:
            dic[site]=[fc]
        else:
            dic[site].append(fc)

    # If the same phosphosite has been detected more than once, a mean log2(FC) is calculated for it.
    # Final array contains unique phosphosites and individual log2(FCs), averaged where appropriate.
    for key in dic:
        length=len(dic[key])
        mean_fc=sum(dic[key])/length
        dic[key] = float(mean_fc)

    # The mean log2(FC) and the standard deviation of all phosphosites are calculated here.
    # These values will be used to obtain the z-score for each kinase later on.
    all_log2=[]    
    for key in dic:
        all_log2.append(dic[key])

    all_mean=sum(all_log2) / float(len(all_log2))
    all_std=np.std(all_log2)

    # Each phosphosite in dic is scanned against the K-S db. 
    # If a match is found, relevant information for that phosphosite is appended to a new array ks_links.
    ks_links=[]
    for x in dic:
        for y in ks_db:
            if x == y[0]:
                ks_links.append([y[1], y[0], y[2], y[3], dic[x]])
            else:
                continue

    # The array is then converted into a dataframe.
    ks_links_df = pd.DataFrame(ks_links, columns=["Kinase", "Site", "Site.Seq(+/- 7AA)", "Source", "log2(FC)"])

    # A dictionary containing unique kinases and substrate log2(FCs) as values is created.
    # If the same kinase was identified for multiple substrates, multiple log2(FCs) are appended to the dictionary values.
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
    # The algorithm then computes the z-score.
    # A new array stores kinase gene, no. of substrates, mean log2(FC), z-score and p-value.
    zscore_info=[]
    for key in kinase_dic:
        substrate_num = len(kinase_dic[key])
        kin_fc_mean = sum(kinase_dic[key]) / float(substrate_num)
        z_score = zScore(kin_fc_mean, all_mean, substrate_num, all_std)
        z_pval = getpValue(z_score)
        zscore_info.append([key, substrate_num, kin_fc_mean, z_score, z_pval])

    # z-score array is converted into a pandas dataframe.
    zscore_df=pd.DataFrame(zscore_info, columns=["Kinase", "Sub.Count", "mnlog2(FC)", "zSc", "pVal"])

    # Barplot is only generated if the user chose to produce graphics at file upload.
    if graphics == "no":
        svg_fig = "Barplot was not generated for this analysis."
    elif graphics == "yes":
        # For easier plot manipulation, the array is sorted in an ascending order according to the z-score.
        zscore_info.sort(key=lambda x: float(x[3]))

        # Kinase label for y-axis.
        kinases = []

        # Corresponding z-scores.
        score_list = []

        # Extraction of desirable results that match user-defined parameters.
        # min_sub denotes the minimum substrate count per kinase indicated by the user.
        for entry in zscore_info:
            kin=entry[0]
            sub=entry[1]
            score = entry[3]
            if sub >= min_sub:
                kinases.append(kin)
                score_list.append(score)

        # Dictionary of kinases as keys and p-values as values. 
        # It will serve as a map for the colour coding of the barplot bars.
        kin_map={}
        for entry in zscore_info:
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
        ax.barh(y_pos, score_list, color=colours)

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
        plt.xlabel("Z-score", fontsize=6)

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
    zscore_df = zscore_df.to_json(orient='split')
    ks_links_df = ks_links_df.to_json(orient='split')
    
    return zscore_df, ks_links_df, svg_fig