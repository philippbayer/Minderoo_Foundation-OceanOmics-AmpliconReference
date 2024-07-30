import os
from Bio import SeqIO
from collections import defaultdict, Counter
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-i', '--infile', help='Name of input fasta downloaded',
        required = True)
args = parser.parse_args()

if __name__ == '__main__':
    all_seqs = []
    removed_classes = defaultdict(int)
    removed = {} # key: sequence id. value: reason why removed.
    vague = {} # key: sequence id. value: vague reason. 
    # need to check whether we have other species for this genus - 
    # if we have, remove vague sp. 
    # if not, keep as Genus-level so we have at least some level of representation
    genus_hits = {}

    long_to_short = {} # key: 12S_NC_12321.1:12321, value: NC_12321.1
    
    if not os.path.exists('0-taxoncheck'):
        os.mkdir('0-taxoncheck')
    with open('0-taxoncheck/TAXONKIT_nuccore.txt', 'w') as out,\
            open('0-taxoncheck/TAXONKIT_gene.txt', 'w') as out2:
        for s in SeqIO.parse(args.infile, 'fasta'):
            # get the taxonomy ID
            # we use entrez
            if ':' in s.id:
                # this is a gene ID that I mangled. Unmangle
                thisid = '_'.join(s.id.split(':')[0].split('_')[1:])
                long_to_short[s.id] = thisid
                out2.write(f'{thisid}\n')
            else:
                out.write(f'{s.id}\n')
            all_seqs.append(s)

    os.popen(f'cat {out.name} | efetch -db nuccore -format docsum > 0-taxoncheck/efetch_nuccore.txt').read()
    os.popen(f'cat 0-taxoncheck/efetch_nuccore.txt | xtract -pattern DocumentSummary -element AccessionVersion,TaxId,Organism > 0-taxoncheck/efetch_names_taxids.txt').read()

    os.popen(f'cat {out2.name} | efetch -db nuccore -format docsum > 0-taxoncheck/efetch_gene.txt').read()
    os.popen(f'cat 0-taxoncheck/efetch_gene.txt | xtract -pattern DocumentSummary -element AccessionVersion,TaxId,Organism >> 0-taxoncheck/efetch_names_taxids.txt').read()

    taxid_dict = {}

    # we also have to check whether we have other species-level
    # vouchers for those entries which are sp./cf.
    # if we have a species-level voucher in a diferent sequence, we can remove the sp.
    genus_dict = defaultdict(set)
    with open('0-taxoncheck/efetch_names_taxids.txt') as fh:
        for line in fh:
            ll = line.rstrip().split('\t')
            try:
                name, taxid, species = ll[0], ll[1], ll[2]
                genus = species.split(' ')[0]
            except:
                continue
            taxid_dict[name] = (taxid, species)
            genus_dict[genus].add(species)

    # now self-blast - we need to make the blast database
    if not os.path.exists('1-selfblast'):
        os.mkdir('1-selfblast')

    with open('1-selfblast/selfblastdb.fasta', 'w') as fastaout\
            , open('1-selfblast/selfblastdb_taxIDs.txt','w') as taxout:
        for s in all_seqs:
            # naming scheme:
            # ID, taxID, speciesname, OLD description

            if ':' in s.id:
                # this is a gene ID that I mangled. Unmangle
                thisid = '_'.join(s.id.split(':')[0].split('_')[1:])
            else:
                thisid = s.id

            taxid, species = taxid_dict[thisid]
            genus = species.split(' ')[0]
            # let's check -is this a vague species?

            # remove vague species
            if 'cf.' in s.description or 'sp.' in s.description or ' x ' in s.description or ' X ' in s.description \
                    or 'aff.' in s.description:
                # get the other species we have for this genus
                all_species_for_genus = genus_dict[genus]
                # we also track all other sp and cf - we have to ignore these
                all_species_for_genus_clean = set()
                for l in all_species_for_genus:
                    if 'cf.' in l or 'sp.' in l or ' x ' in l or ' X ' in l \
                            or 'aff.' in l:
                        continue
                    all_species_for_genus_clean.add(l)

                if len(all_species_for_genus_clean) >= 1: 
                    # we have other species
                    removed[s.id] =  f'vague species identification ({s.description}), {len(all_species_for_genus_clean)-1} other species exist'
                    continue
                # there are some weird cases where people didn't even put the .sp into Genus level, those are sometimes on the
                # family level
                if genus.endswith('dae'):
                    removed[s.id] =  f'vague species identification ({s.description}), with family-level - too high so removing'
                    continue

            new_name = f'{s.id} {taxid} {species} {s.description}'
            fastaout.write(f'>{new_name}\n{str(s.seq)}\n')
            taxout.write(f'{s.id} {taxid}\n')

    # make the blast database
    os.popen('makeblastdb -dbtype nucl -in 1-selfblast/selfblastdb.fasta -parse_seqids -taxid_map 1-selfblast/selfblastdb_taxIDs.txt').read()
    if not os.path.exists('1-selfblast/taxdb.tar.gz'):
        os.popen('wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz -O 1-selfblast/taxdb.tar.gz').read()
        os.popen('tar xzvf 1-selfblast/taxdb.tar.gz --directory 1-selfblast/').read()

    # run blast
    if not os.path.exists('1-selfblast/selfblastdb.results.tsv'):
        os.popen('blastn -db 1-selfblast/selfblastdb.fasta -query 1-selfblast/selfblastdb.fasta -out 1-selfblast/selfblastdb.results.tsv -num_threads 100 -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -perc_identity 100 -qcov_hsp_perc 98').read()

    # now calculate LCA with my script
    os.popen('python computeLCA.py 1-selfblast/selfblastdb.results.tsv > 1-selfblast/selfblastdb.LCAs.tsv').read()
    if not os.path.exists('2-LCAs/'):
        os.mkdir('2-LCAs/')

    # now find the weirdos

    weird_queries = set()

    with open('1-selfblast/selfblastdb.LCAs.tsv') as fh:
        for line in fh:
            ll = line.rstrip().split('\t')
            lineage = ll[1]
            query = ll[0]
            # "{k};{p};{c};{o};{f};{g};{s}"
            # cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Chondrichthyes;Elasmobranchii;Selachii;Galeomorphii;Galeoidea;Carcharhiniformes;Triakidae;Mustelus;Mustelus manazo 
            lineage_split = lineage.split(';')
            # sometimes, there's an empty line
            if lineage_split == ['', '0']: continue

            family, genus, species = lineage_split[-3:]
            if not family and not genus:
                weird_queries.add(query)

    queries_to_weird_subjects = defaultdict(set) # key: weird query, subject: set of all hits for this query
    subjects_to_species = {}
    weird_subjects_to_queries = defaultdict(set)

    with open('1-selfblast/selfblastdb.results.tsv') as fh:
        for line in fh:
            ll = line.split('\t')
            q = ll[0]
            s = ll[1]
            taxid, spec = ll[2], ll[3]
            if q in weird_queries:
                queries_to_weird_subjects[q].add(s)
                weird_subjects_to_queries[s].add(q)
                subjects_to_species[s] = (taxid, spec)

    with open('2-LCAs/TEMP_TAX.txt', 'w') as out:
        for seq_id in subjects_to_species:
            taxid, spec = subjects_to_species[seq_id]
            out.write(f'{seq_id}\t{taxid}\t{spec}\n')

    seq_to_fam = {}

    os.popen('taxonkit lineage -i 2 2-LCAs/TEMP_TAX.txt | taxonkit reformat -i 4 > 2-LCAs/TEMP_TAX_LINEAGE.txt').read()
    with open('2-LCAs/TEMP_TAX_LINEAGE.txt') as fh:
        for line in fh:
            #gb|KC136576.1|  334875  Hyperoglyphe antarctica cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Teleostomi;Euteleostomi;Actinopterygii;Actinopteri;Neopterygii;Teleostei;Osteoglossocephalai;Clupeocephala;Euteleosteomorpha;Neoteleostei;Eurypterygia;Ctenosquamata;Acanthomorphata;Euacanthomorphacea;Percomorphaceae;Pelagiaria;Scombriformes;Centrolophidae;Hyperoglyphe;Hyperoglyphe antarctica Eukaryota;Chordata;Actinopteri;Scombriformes;Centrolophidae;Hyperoglyphe;Hyperoglyphe antarctica
            ll = line.rstrip().split('\t')
            seq_id, tax_id, spec_id, lineage, reform_lineage = ll
            # "{k};{p};{c};{o};{f};{g};{s}"
            kingdom, phylum, thisclass, order, family, genus, species = reform_lineage.split(';')
            seq_to_fam[seq_id] = family
        
    with open('2-LCAs/TEMP_AGAIN.txt', 'w') as out:
        for line in open('1-selfblast/selfblastdb.fasta'):
            if line.startswith('>'):
                ll = line.lstrip('>').split(' ')
                thisid, thistaxid, thisgenus, thisspec = ll[0], ll[1], ll[2], ll[3]
                out.write(f'{thisid}\t{thisgenus} {thisspec}\n')

    os.popen('taxonkit name2taxid -i 2 2-LCAs/TEMP_AGAIN.txt | taxonkit lineage -i 3 > 2-LCAs/TEMP_AGAIN_LINEAGE.txt').read()
    query_to_fam = {}
    for line in open('2-LCAs/TEMP_AGAIN_LINEAGE.txt'):
        ll = line.rstrip().split('\t')
        if len(ll) == 2:
            continue
        seq_id, species_name,taxid,lineage = ll
        # "{k};{p};{c};{o};{f};{g};{s}"
        family = 'NA'
        for element in lineage.split(';'):
            if element.endswith('dae'):
                family = element
        query_to_fam[seq_id] = family

    family_count = Counter()
    seen_iffys = set()
    with open('2-LCAs/Iffy_sequences.tsv', 'w') as out:
        out.write('Sequence ID\tListed family\tShould-be family\n')
        for query in queries_to_weird_subjects:
            all_q = queries_to_weird_subjects[query]
            try:
                x = query_to_fam[query]
            except:
                continue
            family_count = Counter()
            fam_to_seq = defaultdict(list)
            for q in all_q:
                family_count[seq_to_fam[q]] += 1
                fam_to_seq[seq_to_fam[q]].append( q)
            most_common = family_count.most_common(100)
            iffy_ones = most_common[1:]
            majority = most_common[0]
            # [('Scombridae', 5), ('Balistidae', 1)]
            # becomes Balistidae alone
            # which sequence is this?
            iffy_sequences = []
            for i in iffy_ones:
                fam = i[0]
                iffy_sequences += fam_to_seq[fam]
            #print(iffy_sequences, iffy_ones, 'should be ', majority)
            # ['gb|HQ592314.1|', 'gb|HQ592315.1|', 'gb|HQ592313.1|'] [('Ophidiidae', 3)] should be  ('Centrolophidae', 5)
            for i in iffy_sequences:
                if i not in seen_iffys:
                    if '|' in i:
                        shorti = i.split('|')[1]
                    else:
                        shorti = i
                    removed[shorti] = f'Placed as {seq_to_fam[i]}, might be {majority[0]}'
                    seen_iffys.add(i)

    if not os.path.exists('3-Final'):
        os.mkdir('3-Final')

    with open('3-Final/Final_database.fasta', 'w') as outfasta, \
            open('3-Final/QC_Removal_stats.tsv', 'w') as out,\
            open('3-Final/Final_database_taxids.txt', 'w') as outtaxa:
        for i in removed:
            out.write(f'{i}\t{removed[i]}\n')
        for seq in SeqIO.parse('1-selfblast/selfblastdb.fasta', 'fasta'):
            if seq.id in removed:
                continue
            thisid, thistaxid, *bla = seq.description.split(' ')
            outtaxa.write(f'{thisid} {thistaxid}\n')
            outfasta.write(seq.format('fasta'))
    
    # now format the BLAST database
    os.popen('makeblastdb -dbtype nucl -in 3-Final/Final_database.fasta -parse_seqids -taxid_map 3-Final/Final_database_taxids.txt').read()

    # now make stats for removal
    for i in removed:
        removed_classes[removed[i].split('(')[0]] += 1

    print('Final output is in "3-Final/Final_database.fasta".')
    print('Final stats are in "3-Final/QC_Removal_stats.tsv".')
    print('Summary stats are in "3-Final/QC_Removal_summary_stats.tsv".')
    print(f'Overall stats: removed {sum(removed_classes.values())} sequences.')
    with open("3-Final/QC_Removal_summary_stats.tsv", 'w') as out:
        for c in removed_classes:

            out.write(f'{c}\t{removed_classes[c]}\n')
