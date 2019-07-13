# A parsed version of the PSP database is uploaded.
def uploadPSP():
    ks_data=[]
    pp_db=open("databases/psp_db.tsv", "r")
    headers=pp_db.readline()
    for line in pp_db:
        line=line.rstrip("\n")
        line=line.split("\t")
        sub=line[0]
        kin=line[3]
        site_seq=line[6]
        psp_source=line[7]
        ks_data.append((sub, kin, site_seq, psp_source))
    pp_db.close()
    
    # In case the PSP database contains duplicates, these are removed here.
    # Duplicate phosphosite entries are removed even if the reported sequence is different for each. Only one copy is retained.
    #unique_ks=set(ks_data)
    #ks_db=list(unique_ks)
    unique_ks=dict(((x[0], x[1]), x) for x in ks_data).values()
    psp_db=list(unique_ks)
    
    return psp_db

# PDTS database is uploaded. 
# Duplicates have already been pre-removed in the script generating the .tsv based database.
def uploadPDTS():
    pdtsdb=open("databases/pdts_db.tsv", "r")
    header=pdtsdb.readline()
    pdts_db=[]
    for line in pdtsdb:
        line=line.rstrip("\n")
        line=line.split("\t")
        kin=line[0]
        sub=line[1]
        source=line[5]
        pdts_db.append([sub, kin, source])
    pdtsdb.close()
    
    return pdts_db

# EDGES database is uploaded.
# Duplicates have been removed in the script generating the database file.
def uploadEDGES():
    edges = open("databases/edges_db.tsv", "r")
    header = edges.readline()
    edges_db = []
    for line in edges:
        line = line.rstrip("\n")
        line = line.split("\t")
        site = line[0]
        kin = line[1]
        source = line[2]
        edges_db.append([site, kin, source])
    edges.close()
    
    return edges_db