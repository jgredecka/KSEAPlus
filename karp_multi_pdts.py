def userInput(df, min_sub, ks_db, graphics):

    import io
    from io import StringIO
    import pandas as pd
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

    # Function to calculate total number of known phosphosites in the PDTS database for a given kinase.
    def getTotalSub(kinase):
        counter=0
        for x in ks_db:
            db_kinase = x[1]
            if kinase == db_kinase:
                counter+=1
        return counter

    # Function to calculate the K-Score.
    def kScore(kin_sum, all_sum, sub_num, total_sub):
        kscore = (kin_sum/all_sum) * (sub_num/total_sub)**(1/2) * 10**6
        return kscore

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

    dic={}
    kinase_dic={}
    kscore_info=[]
    ks_links=[]
    ks_info=[]
    heatmap_array=[]
    kscore_columns=["Kinase", "Sub.Count", "Total.Sub.Count"]
    ks_columns=["Kinase", "Site", "Source"]
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
        kscore_columns.append("Sum.Ints." + curr_col)
        kscore_columns.append("kSc." + curr_col)
        ks_columns.append("Ints." + curr_col)
        heatmap_col.append(curr_col)

        data=[]
        all_ints=[]
        all_ints=[]
        # Multiple phosphosites separated by a colon are split here.
        # This ensures each phosphosite substrate and the intensity value starts with a new line.
        # This is ran for each sample in turn.
        for n in range (0, len(array)):
            site=array[n][0].upper()
            ints=array[n][col]
            site=site.split(";")
            for s in site:
                if s == '':
                    continue
                else:
                    # If "no_mod" is found in any genes, the phosphosite is ambiguous and is omitted.
                    if "NO_MOD" not in s:
                        data.append([s, float(ints)])

        # Mapping of phosphosite substrate keys to their (often multiple) intensity values is achieved here.       
        for entry in data:
            site=entry[0]
            ints_val=entry[1]
            if site not in dic:
                dic[site]=[ints_val]
            else:
                dic[site].append(ints_val)

        # If the same phosphosite has been detected more than once, its mean inetensity is calculated.
        # Final dictionary contains unique phosphosites and individual intensity values, averaged where appropriate.
        for key in dic:
            length=len(dic[key])
            mean_ints=sum(dic[key])/length
            dic[key] = float(mean_ints)

        # For each sample, the sum of all phosphosite intensities in the dataset is calculated here. This will be used to obtain a K-score for each identified kinase later on.
        for key in dic:
            all_ints.append(dic[key])
        sum_ints=sum(all_ints)

        # Each phosphosite in the dictionary is scanned against the PDTS K-S db. 
        # If a match is found, relevant information for that phosphosite is retained.
        # Scanning is only done for the first column.
        if col == 1:
            for x in dic:
                for y in ks_db:
                    if x == y[0]:
                        # ks_links will be used to assign the current sample's intensity to each kinase later on.
                        ks_links.append([y[1], y[0], y[2], dic[x]])
                        # ks_info will contain kinase-substrate relationship info for each sample.
                        ks_info.append([y[1], y[0], y[2], dic[x]])
        # Once the first column is passed, new intensities are removed and/or appended to the original arrays for each sample.
        elif col > 1:
            for s in ks_links:
                s.remove(s[-1])
                s.append(dic[s[1]])
            for k in ks_info:
                k.append(dic[k[1]])

        # A dictionary containing unique kinases and substrate intensities is created.
        # If the same kinase was identified for multiple substrates, multiple intensities are appended to the dictionary values.
        for match in ks_links:
            kinase=match[0]
            intensity=match[3]
            if kinase not in kinase_dic:
                kinase_dic[kinase]=[intensity]
            else:
                kinase_dic[kinase].append(intensity)

        # The dictionary is used to calculate the number of substrates identified for each unique kinase.
        # It also calculates the sum of intensities across each kinase's substrates.
        # A predefined function computes the total substrate count in the DB for each kinase.
        # The algorithm then computes the k-score.
        # All information is appended to a new array for each sample.
        condition_ind=-1
        index=-1
        for kinase in kinase_dic:
            index+=1
            # Number of phosphosites in the dataset identified for a given kinase.
            sub_num=len(kinase_dic[kinase])
            # Sum of intensities of phosphosites associated with a given kinase.
            kin_ints_sum = sum(kinase_dic[kinase])
            # Total number of sites in the DB for a given kinase.
            total_sub = getTotalSub(kinase)
            # Kinase 'K-Score' is calculated here.
            kscore = kScore(kin_ints_sum, sum_ints, sub_num, total_sub)
            if col == 1:
                # KSEA results for each kinase. 
                kscore_info.append([kinase, sub_num, total_sub, kin_ints_sum, kscore])
            else:
                # If the program has gone past the first sample, kin_ints_sum and k_score are appended in a repeating manner to the original array.
                kscore_info[index].append(kin_ints_sum)
                kscore_info[index].append(kscore)
            # For the heatmap, only k-scores for the kinases with the minimum substrate count specified by the user are extracted.
            if sub_num >= min_sub:
                condition_ind+=1
                if col == 1:
                    # An array of kinases and k-scores corresponding to multiple samples. Used for the heatmap.
                    heatmap_array.append([kinase, kscore])
                else:
                    heatmap_array[condition_ind].append(kscore)


    # Score and kinase-substrate relationships DFs are generated.
    kscore_df = pd.DataFrame(kscore_info, columns=kscore_columns)
    ks_df = pd.DataFrame(ks_info, columns=ks_columns)
    # Heatmap df for the heatmap generation.
    heatmap_df = pd.DataFrame(heatmap_array, columns=heatmap_col)
    heatmap_df=heatmap_df.set_index("Kinase")

    # Heatmap only generated if the user chose to produce graphics during file upload.
    if graphics == "no":
        svg_fig = "Heatmap was not generated for this analysis."
    elif graphics == "yes":
        # Set the margins and square height for a single category.
        topmargin = 0.1 #inches
        bottommargin = 0.1 #inches
        categorysize = 0.25 # inches
        # Number of kinases identified.
        n=len(heatmap_array)
        leftmargin = 0.1
        rightmargin = 0.1
        catsize = 0.3
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
        ax = sns.heatmap(heatmap_df, cmap='coolwarm', cbar=False, linewidths=0.3, linecolor='white')

        # Format the colour bar dynamically.
        ax_div = make_axes_locatable(ax)
        width = axes_size.AxesY(ax, aspect=1./aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax = ax_div.append_axes('right', size = width, pad = pad)
        cb=plt.colorbar(ax.get_children()[0], cax = cax, orientation = 'vertical')
        cax.yaxis.set_ticks_position('right')
        cb.ax.tick_params(labelsize=6)
        cb.set_label('K-Score', fontsize=6, labelpad=7)
        cb.outline.set_visible(False)

        #Remove y axis label.
        ax.yaxis.set_label_text("")

        # Rotate the axis labels.
        for item in ax.get_yticklabels():
            item.set_rotation(0)

        for item in ax.get_xticklabels():
            item.set_rotation(90)

        # Create a StringIO object and use it to write SVG figure data to string buffer.
        fig_file = StringIO()
        fig.savefig(fig_file, format='svg', bbox_inches="tight")
        # Seek beginning of the figure file.
        fig_file.seek(0)
        # Retrieve figure contents as a string.
        svg_fig = '<svg' + fig_file.getvalue().split('<svg')[1]
        # Free memory buffer.
        fig_file.close()
    
    return kscore_df, ks_df, svg_fig