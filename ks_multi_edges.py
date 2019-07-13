def userInput(df, min_sub, ks_db, graphics):

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

    # Function that returns the Kolmogorov-Smirnov test statistic along with the p-value.
    # -log10 of p-value is also calculated for the barplot.
    def ksTest(substrates, non_substrates):
        result = st.ks_2samp(substrates, non_substrates)
        ks_stat = result[0]
        pval = result[1]
        log_pval = np.log10(1/pval)
        return ks_stat, pval, log_pval

    # User data is passed from the server and parsed as appropriate.
    user_file = df.values.tolist()
    header = df.columns.values.tolist()
    col_length=len(header)
    array=[]
    for line in user_file:
        if line == "":
            continue
        else:       
            if "" not in line:
                array.append(line)

    nonsub_dic={}
    kinase_dic={}
    dic={}
    pval_map=[]
    heatmap_array=[]
    kolsmir_info=[]
    ks_links=[]
    ks_info=[]
    kolsmir_col=["Kinase", "Sub.Count"]
    ks_col=["Kinase", "Site", "Source"]
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
        kolsmir_col.append("mnlog2(FC)." + curr_col)
        kolsmir_col.append("(+/-)KS." + curr_col)
        kolsmir_col.append("pVal." + curr_col)
        kolsmir_col.append("(+/-)-log10(pVal)." + curr_col)
        ks_col.append("log2(FC)." + curr_col)
        heatmap_col.append(curr_col)

        data=[]
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
            
        # Each phosphosite in the dictionary is scanned against the EDGES K-S db. 
        # If a match is found, relevant information for that phosphosite is retained.
        # Scanning is only done for the first column.
        if col == 1:
            for x in dic:
                for y in ks_db:
                    if x == y[0]:
                        # ks_links will be used to assign the current sample's log2(FCs) to each kinase later on.
                        ks_links.append([y[1], y[0], y[2], dic[x]])
                        # ks_info will contain kinase-substrate relationship info for each sample.
                        ks_info.append([y[1], y[0], y[2], dic[x]])
        # Once the first column is passed, new log2(FCs) are removed and/or appended to the original arrays for each sample.
        elif col > 1:
            for s in ks_links:
                s.remove(s[-1])
                s.append(dic[s[1]])
            for k in ks_info:
                k.append(dic[k[1]])
        
        # List converted into a dataframe for further data manipulation.            
        ks_links_df = pd.DataFrame(ks_links, columns = ["Kinase", "Site", "Source", "log2(FC)"])
        
        # A dictionary containing unique kinases and substrate log2(FCs) is created.
        # If the same kinase was identified for multiple substrates, multiple log2(FCs) are appended to the dictionary values.
        for match in ks_links:
            kinase=match[0]
            log2fc=match[3]
            if kinase not in kinase_dic:
                kinase_dic[kinase]=[log2fc]
            else:
                kinase_dic[kinase].append(log2fc)
                
        # The dictionary is used to calculate the number of substrates identified for each unique kinase.
        # It also calculates the mean log2(FC) across each kinase's substrates.
        # The algorithm computes the (+/-)KS statistic, p-value and -log10(p-value) using substrate and non-substrate log2(FC) values.
        index=-1
        condition=-1
        for kinase in kinase_dic:
            index+=1
            nonsub_fc=[]
            sub_num=len(kinase_dic[kinase])
            kin_fc_mean=sum(kinase_dic[kinase]) / float(sub_num)
            if col == 1:
            # Substrate names for a given kinase are extracted into a flat list.
            # These will be used to identify all non-substrates conditionally.
                sub_df = ks_links_df.loc[ks_links_df['Kinase'] == kinase]
                sub_names = sub_df['Site'].tolist()
                # Non-substrates are located within the ks_links dataframe.
                # For faster performance, these are assigned to a new dictionary to be re-used for columns 2 onwards.
                non_sub_df=ks_links_df.loc[~ks_links_df['Site'].isin(sub_names)]
                non_sub_df=non_sub_df.drop_duplicates(subset='Site')
                nonsub_names = non_sub_df['Site'].tolist()
                nonsub_dic[kinase] = nonsub_names
                # Substrate and non-substrate log2(FC) values.
                sub_fc = kinase_dic[kinase]
                non_sub_fc = non_sub_df['log2(FC)'].tolist()
                #KS-test statistic and p-value for each kinase are calculated here.
                ks_stat, pval, log_pval = ksTest(sub_fc, non_sub_fc)
                # -log10(p-val) and KS is signed based on the mean log2(FC) of the kinase.
                if kin_fc_mean < 0:
                    log_pval = -log_pval
                    ks_stat = -ks_stat
                # Kolmogorov-Smirnov stats are appended to a new list here.
                kolsmir_info.append([kinase, sub_num, kin_fc_mean, ks_stat, pval, log_pval])

            elif col > 1:
                # Non-substrate log2(FCs) are identified for each kinase.
                sub_fc = kinase_dic[kinase]
                non_subs = nonsub_dic[kinase]
                for n in non_subs:
                    nonsub_fc.append(dic[n])
                # KS-test function applied here.
                ks_stat, pval, log_pval = ksTest(sub_fc, nonsub_fc)
                # -log10(p-val) and KS is signed based on the mean log2(FC) of the kinase.
                if kin_fc_mean < 0:
                    log_pval = -log_pval
                    ks_stat = -ks_stat
                # If the program has gone past the first column, each statistic is appended in a repeating manner to the original array.
                kolsmir_info[index].append(kin_fc_mean)
                kolsmir_info[index].append(ks_stat)
                kolsmir_info[index].append(pval)
                kolsmir_info[index].append(log_pval)

            if sub_num >= min_sub:
                condition+=1
                if col == 1:
                    # Array used for heatmap generation.  
                    heatmap_array.append([kinase, log_pval])
                    # An array of p-values for each kinase across all samples. Used for heatmap annotation.
                    pval_map.append([pval])
                elif col > 1:
                    heatmap_array[condition].append(log_pval)
                    pval_map[condition].append(pval)

    # p-values for the heatmap annotation are extracted from a nested list into a flast list.
    pvalues=[]
    for entry in pval_map:
        for pval in entry:
            pvalues.append(pval)

    # Statistic and KS-links dataframes are generated.
    kolsmir_df = pd.DataFrame(kolsmir_info, columns=kolsmir_col)
    ksinfo_df = pd.DataFrame(ks_info, columns=ks_col)
    # Heatmap df for the heatmap generation.
    heatmap_df = pd.DataFrame(heatmap_array, columns=heatmap_col)
    heatmap_df=heatmap_df.set_index("Kinase")

    # Heatmap only generated if the user chose to produce graphics during file upload.
    if graphics == "no":
        svg_fig = "Heatmap was not generated for this analysis."
    elif graphics == "yes":
        # Set the margins and bar height for a single category.
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

        # Calculate a dynamic figure height based on the known values above.
        figheight = topmargin + bottommargin + (n+1)*categorysize

        # Calculate a dynamic figure width based on the known values above.
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
        cb.set_label('(+/-) -log10(p-value)', fontsize=6, labelpad=7)
        cb.outline.set_visible(False)

        #Remove y axis label.
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
    
    return kolsmir_df, ksinfo_df, svg_fig
    