def userInput(df, min_sub, ks_db, graphics):

    import io
    from io import StringIO
    import pandas as pd
    import numpy as np
    
    if graphics == "yes":
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt

    # Function to calculate total number of known phosphosites in the database for a given kinase.
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
    user_file=df.values.tolist()
    array=[]
    for line in user_file:
        phosphosite=line[0]
        intensity = line[1]
        if phosphosite == "" or intensity == "":
            continue
        else:
            array.append([phosphosite, float(intensity)])

    # Multiple phosphosites separated by a colon are split here.
    # This ensures each phosphosite substrate starts with a new line.
    data=[]
    for n in range (0, len(array)):
        site=array[n][0].upper()
        ints=array[n][1]
        site=site.split(";")
        for s in site:
            if s == '':
                continue
            else:
                # If "no_mod" is found in any genes, the phosphosite is ambiguous and is omitted.
                if "NO_MOD" not in s:
                    data.append([s, ints])

    # Mapping of phosphosites to their (often multiple) intensity values is achieved here.
    dic={}
    for entry in data:
        site=entry[0]
        ints_val=entry[1]
        if site not in dic:
            dic[site]=[ints_val]
        else:
            dic[site].append(ints_val)

    # If the same phosphosite has been detected more than once, a mean intensity is calculated for that phosphosite.
    # Final array contains unique phosphosites and individual intensity values, averaged where appropriate.
    for key in dic:
        length=len(dic[key])
        mean_ints=sum(dic[key])/length
        dic[key] = float(mean_ints)

    # The sum of all phosphosite intensities in the dataset is calculated here.
    # This will be used to calculate the K-score for each kinase later on.
    all_ints=[]
    for key in dic:
        all_ints.append(dic[key])
    sum_ints=sum(all_ints)

    # Each phosphosite substrate in dic is scanned against the PSP K-S db. 
    # If a match is found, relevant information for that phosphosite is appended to a new array called ks_links.
    ks_links=[]
    for x in dic:
        for y in ks_db:
            if x == y[0]:
                ks_links.append([y[1], y[0], y[2], y[3], dic[x]])
            else:
                continue

    # The array is then converted into a dataframe to be viewed as a table.
    ks_links_df = pd.DataFrame(ks_links, columns=["Kinase", "Site", "Site.Seq(+/- 7AA)", "Source", "Ints"])

    # A dictionary containing unique kinases as keys and substrate intensities as values is created.
    # If the same kinase was identified for multiple substrates, multiple intensities are appended to the dictionary values.
    kinase_dic={}
    for match in ks_links:
        kinase=match[0]
        intensity=match[4]
        if kinase not in kinase_dic:
            kinase_dic[kinase]=[intensity]
        else:
            kinase_dic[kinase].append(intensity)

    # The dictionary is used to calculate the number of substrates identified for each kinase.
    # It also calculates the sum of intensities across each kinase's substrates.
    # A predefined function computes the total substrate count in the PSP db for each kinase.
    # The algorithm then computes the k-score.
    # All information is appended to a new array.
    kscore_info=[]
    for kinase in kinase_dic:
        # Number of phosphosites in the dataset identified for a given kinase.
        substrate_num = len(kinase_dic[kinase])
        # Sum of intensities of phosphosites associated with a given kinase.
        kin_ints_sum = sum(kinase_dic[kinase])
        # Total number of sites in the PSP database for a given kinase.
        total_sub_num = getTotalSub(kinase)
        # Kinase 'K-Score' is calculated here.
        kscore = kScore(kin_ints_sum, sum_ints, substrate_num, total_sub_num)
        # An array contains all relevant k-score information for each kinase. 
        kscore_info.append([kinase, substrate_num, total_sub_num, kin_ints_sum, kscore])

    # k-score information array is converted into a dataframe.
    kscore_df=pd.DataFrame(kscore_info, columns=["Kinase", "Sub.Count", "Total.Sub.Count", "Sum.Ints", "kSc"])

    # Barplot only generated if the user chose to produce graphics during file upload.
    if graphics == "no":
        svg_fig = "Barplot was not generated for this analysis."
    elif graphics == "yes":
        # The array is sorted in an ascending order to make list manipulaton for plot generation easier.
        kscore_info.sort(key=lambda x: float(x[4]))

        # Kinase label for y-axis.
        kinases = []

        # Corresponding k-scores.
        score_list = []

        # Extraction of desirable results that match user_defined parameters.
        # min_sub denotes the minimum substrate count per kinase indicated by the user.
        for entry in kscore_info:
            kin=entry[0]
            sub=entry[1]
            score = entry[4]
            if sub >= min_sub:
                kinases.append(kin)
                score_list.append(score)

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
        ax.barh(y_pos, score_list, color="#b3b3b3")

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
        plt.xlabel("K-score", fontsize=6)

        # Create a StringIO object and use it to write SVG figure data to string buffer.
        fig_file = StringIO()
        fig.savefig(fig_file, format='svg', bbox_inches="tight")
        # Seek beginning of the figure file.
        fig_file.seek(0)
        # Retrieve figure contents as a string.
        svg_fig = '<svg' + fig_file.getvalue().split('<svg')[1]
        # Free memory buffer.
        fig_file.close()
    
    return kscore_df, ks_links_df, svg_fig