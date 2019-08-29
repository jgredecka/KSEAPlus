def userInput(ks_db, graphics, df, min_sub):
    
    import io
    from io import StringIO
    import pandas as pd
    import scipy.stats as st
    import numpy as np
    
    if graphics == "yes":
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import axes_size
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
        from mpl_toolkits.axes_grid1.axes_size import AxesY, Fraction
        from mpl_toolkits.axes_grid1.colorbar import colorbar
        import seaborn as sns
    
    # Function used to convert kinase z-scores to corresponding p-values.
    # It is assmued z-scores are normally distributed.
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
    f = pd.read_json(df, orient="split")
    
    # User data is passed from the server and parsed as appropriate.
    user_file = f.values.tolist()
    header = f.columns.values.tolist()
    col_length=len(header)
    array=[]
    for line in user_file:
        if line == "":
            continue
        else:       
            if "" not in line:
                array.append(line)
    
    kinase_dic={}
    dic={}
    ks_links=[]
    ks_info=[]
    pval_map=[]
    heatmap_array=[]
    zscore_info=[]
    z_columns=["Kinase", "Sub.Count"]
    ks_columns=["Kinase", "Site", "Site.Seq(+/- 7AA)", "Source"]
    heatmap_col=["Kinase"]
    # Columns 1 and onwards represent samples (e.g. cell lines).
    # For the current column (sample) a set of operations is performed. 
    for col in range(1, col_length):
        # Reset the values in each dictionary for a new column.
        for key in dic:
            dic[key] = []
        for kin in kinase_dic:
            kinase_dic[kin] = []
        # Column names for relevant dataframes are created here dynamically.
        # curr_col is current sample/column name.
        curr_col=header[col]
        z_columns.append("mnlog2(FC)." + curr_col)
        z_columns.append("zSc." + curr_col)
        z_columns.append("pVal." + curr_col)
        ks_columns.append("log2(FC)." + curr_col)
        heatmap_col.append(curr_col)

        data=[]
        all_log2=[]
        # Multiple phosphosites separated by a colon are split here.
        # This ensures each phosphosite substrate and the log2(FC) value starts with a new line.
        # This is ran for each sample in turn. 
        for n in range (0, len(array)):
            site=array[n][0].upper()
            fc=array[n][col]
            site=site.split(";")
            for s in site:
                if s == '':
                    continue
                else:
                    data.append([s, float(fc)])

        # Mapping of phosphosite substrate keys to their (often multiple) log2(FCs) is achieved here.
        for entry in data:
            site=entry[0]
            fc=entry[1]
            if site not in dic:
                dic[site]=[fc]
            else:
                dic[site].append(fc)

        # If the same phosphosite has been detected more than once, its mean log2(FC) is calculated.
        # Final dictionary contains unique phosphosites and individual log2(FC) values, averaged where appropriate.
        for key in dic:
            length=len(dic[key])
            mean_fc=sum(dic[key])/length
            dic[key] = float(mean_fc)

        # For each sample, the mean log2(FC) and standard deviation of all phosphosites in the dataset are calculated here. These values will be used to obtain a z-score for each identified kinase later on.
        for key in dic:
            all_log2.append(dic[key])
        all_mean=sum(all_log2) / float(len(all_log2))
        all_std=np.std(all_log2)

        # Each phosphosite in the dictionary is scanned against the K-S db. 
        # If a match is found, relevant information for that phosphosite is retained.
        # Scanning is only done for the first column.
        if col == 1:
            for x in dic:
                for y in ks_db:
                    if x == y[0]:
                        # ks_links will be used to assign the current sample's log2(FCs) to each kinase later on.
                        ks_links.append([y[1], y[0], y[2], y[3], dic[x]])
                        # ks_info will contain kinase-substrate relationship info for each sample.
                        ks_info.append([y[1], y[0], y[2], y[3], dic[x]])
        # Once the first column is passed, new log2(FCs) are removed and/or appended to the original arrays for each sample.
        elif col > 1:
            for s in ks_links:
                s.remove(s[-1])
                s.append(dic[s[1]])
            for k in ks_info:
                k.append(dic[k[1]])

        # A dictionary containing unique kinases and substrate log2(FCs) is created.
        # If the same kinase was identified for multiple substrates, multiple log2(FCs) are appended to the dictionary values.
        for match in ks_links:
            kinase=match[0]
            log2fc=match[4]
            if kinase not in kinase_dic:
                kinase_dic[kinase]=[log2fc]
            else:
                kinase_dic[kinase].append(log2fc)

        # The dictionary is used to calculate the number of substrates identified for each kinase.
        # It also calculates the mean log2(FC) across each kinase's substrates.
        # The algorithm then computes the z-score.
        # A new array stores kinase gene, no. of substrates, mean log2(FC), z-score and p-value.
        condition_ind=-1
        index=-1
        for key in kinase_dic:
            index+=1
            substrate_num=len(kinase_dic[key])
            kin_fc_mean=sum(kinase_dic[key]) / float(substrate_num)
            z_score = zScore(kin_fc_mean, all_mean, substrate_num, all_std)
            z_pval=getpValue(z_score)
            # For the heatmap, only z-scores for the kinases with the minimum substrate count specified by the user are extracted.
            if substrate_num >= min_sub:
                condition_ind+=1
                if col == 1:
                    # An array of kinases and z-scores corresponding to multiple samples. Used for the heatmap.
                    heatmap_array.append([key, z_score])
                    # An array of p-values for each z-score across all samples. Used for heatmap annotation.
                    pval_map.append([z_pval])
                else:
                    heatmap_array[condition_ind].append(z_score)
                    pval_map[condition_ind].append(z_pval)
            if col == 1:
                # KSEA results stored here.
                zscore_info.append([key, substrate_num, kin_fc_mean, z_score, z_pval])
            # If the program has gone past the first condition column, kin_fc_mean, z_score and zpval are appended in a repeating manner to the original array.
            else: 
                zscore_info[index].append(kin_fc_mean)
                zscore_info[index].append(z_score)
                zscore_info[index].append(z_pval)
                
    # p-values for the heatmap annotation are extracted from a nested list into a flat list.
    pvalues=[]
    for entry in pval_map:
        for pval in entry:
            pvalues.append(pval)            

    # Score and kinase-substrate relationships DFs are generated.
    zscore_df = pd.DataFrame(zscore_info, columns=z_columns)       
    ks_df = pd.DataFrame(ks_info, columns=ks_columns)
    # Heatmap DFs for the heatmap generation.
    heatmap_df = pd.DataFrame(heatmap_array, columns=heatmap_col)
    heatmap_df=heatmap_df.set_index("Kinase")
    
    # Heatmap only generated if the user chose to produce graphics during file upload.
    if graphics == "no":
        svg_fig = "Heatmap was not generated for this analysis."
    elif graphics == "yes":
        # Set the margins and square height for a single category.
        topmargin = 0.1 #inches
        bottommargin = 0.1 #inches
        categorysize = 0.35 # inches
        # Number of kinases identified.
        n=len(heatmap_array)
        leftmargin = 0.1
        rightmargin = 0.1
        catsize = 0.5
        # Number of conditions (e.g. cell lines).
        m=len(heatmap_col)-1

        # Parameters for color bar.
        aspect = n
        pad_fraction = 0.7

        # Calculate a dynamic figure height.
        figheight = topmargin + bottommargin + (n+1)*categorysize

        # Calculate a dynamic figure width.
        figwidth = leftmargin + rightmargin + (m+1)*catsize

        fig, ax = plt.subplots(figsize=(figwidth, figheight))

        # Format the axes.
        ax.xaxis.set_ticks_position('top')
        plt.yticks(fontsize=6)
        plt.xticks(fontsize=6)

        # Plot the heatmap.
        ax = sns.heatmap(heatmap_df, cmap='coolwarm', annot=True, fmt=".1f", annot_kws={'size':5}, cbar=False, linewidths=0.3, linecolor='white')

        # Format the colour bar dynamically.
        ax_div = make_axes_locatable(ax)
        width = axes_size.AxesY(ax, aspect=1./aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax = ax_div.append_axes('right', size = width, pad = pad)
        cb=plt.colorbar(ax.get_children()[0], cax = cax, orientation = 'vertical')
        cax.yaxis.set_ticks_position('right')
        cb.ax.tick_params(labelsize=6)
        cb.set_label('Z-Score', fontsize=6, labelpad=7)
        cb.outline.set_visible(False)

        #Remove y-axis label.
        ax.yaxis.set_label_text("")

        # Rotate the axis labels.
        for item in ax.get_yticklabels():
            item.set_rotation(0)

        for item in ax.get_xticklabels():
            item.set_rotation(90)

        # Annotate statistically significant scores with asterisks.
        # * for p < 0.05 and ** for p < 0.01.
        counter=-1
        for text in ax.texts:
            counter+=1
            if pvalues[counter] < 0.05 and pvalues[counter] >= 0.01:
                text.set_weight('bold')
                text.set_text(text.get_text() + "*")
            elif pvalues[counter] < 0.01:
                text.set_weight('bold')
                text.set_text(text.get_text() + "**")

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
    ks_df = ks_df.to_json(orient='split')
        
    return zscore_df, ks_df, svg_fig