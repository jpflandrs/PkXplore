module Analysis
using Dates
using Random
using SplitApplyCombine
using Luxor
using Base.Iterators
using Serialization

export uniqueutilisateur,faitblast,alignement,trim,faitrapidnj,validate_fasta_nucleic,validate_fasta_proteic,panoramatographe_nuc,exemples,compresser,retourafasta,tranchedemsa




function validate_fasta_nucleic(sequence::String, is_nucleic::Bool)
    nucleic_acid_pattern = r"^>(.+)\r?\n([ATGCNWSMKRYBDHV\-\n]+)$"i
    pattern = is_nucleic ? nucleic_acid_pattern : amino_acid_pattern
    return match(pattern, sequence) !== nothing
end

function validate_fasta_proteic(sequence::String, is_proteic::Bool)
    amino_acid_pattern = r"^>(.+)\r?\n([ACDEFGHIKLMNPQRSTVWYUOBZX\*\-\n]+)$"i
    pattern = is_proteic ? amino_acid_pattern : nucleic_acid_pattern
    return match(pattern, sequence) !== nothing
end

function uniqueutilisateur()
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    random_string = randstring(8)  # 8-character random string
    fichtempo =  joinpath(pwd(),"public","utilisateurs","task_$(timestamp)_$(random_string)")
    mkdir(fichtempo)
    atelier=joinpath(fichtempo,"atelier_"*timestamp*"_"*random_string)
    mkdir(atelier)
    #resu=faitblast(data,fichtempo)

    return atelier
end

function compresser(classeur_utilisateur) #intermédiaire de zippp (oui 3 p) pour travailler dans le directory utilisateur 
    #pwd(),"public","utilisateurs","task_$(timestamp)_$(random_string)"
    
    println("---")
    println(classeur_utilisateur)#/Users/jean-pierreflandrois/Documents/PKDBGENIESTIPPLE/public/utilisateurs/task_20241013_220057_lfXiwPpZ/
    originaldir=pwd()
    println(originaldir)

    cd(classeur_utilisateur)

    
    println(pwd())
    utilisateur=splitdir(classeur_utilisateur)[2]*".tar.zip"
    latelier=replace(splitdir(classeur_utilisateur)[2],"task_" => "atelier_")
    println(utilisateur, "  <- tar va faire <-  ",latelier)
    try 
        cmd=`tar -zcvf  $utilisateur $latelier`
        println(cmd)
        run(pipeline(cmd,stdout=devnull,stderr=devnull))#,stdout="dev/null",stderr="dev/null"))
        
    catch
        println(" ERREUR   targz ",latelier,"   ")
        
    end

    cd(originaldir)
    println("retour...",pwd())
    return utilisateur

end

function faitblast(db_blast::String, S::String, dirutilisateur::String, nbsearch::Int64) #nombnaque sequence dir utilisateur nb items a chercher
    #sauvegarder le fasta d'entrée
    fichiersorties::String=joinpath(dirutilisateur,"blastfasta.fasta")
    write(fichiersorties,S)
    la_db::String = joinpath("BLAST",db_blast) #reconstruire le chemin 
    #blaster
    #la_db="/Users/jean-pierreflandrois/Documents/MVCGenieApp/public/BLAST/TRECS_16SrRNA.fst"
    seqfile::String=joinpath(dirutilisateur,"blastfasta.fasta")
    #db_blast="TRECS_16SrRNA"
    limitedsearch::String= string(nbsearch)

    vecteurtetes::Vector{String}=[]
    vecteurtetescoupees::Vector{String}=[]
    scoredivers::Vector{Tuple{SubString{String}, SubString{String}}}=[]#::Tuple{SubString{String},SubString{String}}=[]
    collection_avec_query::String = ""

    outfile=seqfile*"-blast"*db_blast*"_out.html"
    sorties=joinpath(dirutilisateur,"blast.out.html")
    cmd=`blastn -query $seqfile -db $la_db -evalue 0.1 -dust no -subject_besthit -num_threads 8 -max_target_seqs $limitedsearch -outfmt 0 -html`
    run(pipeline(cmd,stdout=sorties,stderr=devnull)) #stdout="dev/null",stderr="dev/null")
    listedesaextraire=[]
    
    fileis = open(sorties) do f
        acontinuer=true
        while  acontinuer
            for l in eachline(f) #note : if startswith(l,'>') # un petitipetit peu plus long
                if occursin("***** No hits found *****",l)
                    acontinuer=false
                    return vecteurtetescoupees,vecteurtetes,scoredivers,collection_avec_query
                elseif occursin("<a href=#",l)
                    push!(listedesaextraire,split(l,' ')[1])
                    scorecrude=split(l,'>')  
                    #GCF_000430995.1 Thermobrachium_celere~RTES~GCF_000430995.1=Bacter...  <a href=#GCF_000430995.1>1218</a>    0.0
                    push!(scoredivers,(strip(split(scorecrude[2],'<')[1]) ,strip(scorecrude[3]) )) #tuple ("1218","0")
                elseif occursin("<a name=",l)
                    acontinuer=false
                end
            end
        end
    end
    ciblestxt::String=joinpath(dirutilisateur,"blast.list.txt")
    write(ciblestxt,join(listedesaextraire,"\n"))

    f = open(seqfile, "r")
    fastainit::String= read(f, String)       
    close(f)
    extraitsbnk::String=joinpath(dirutilisateur,"collection_de_blast.fasta")
    cmd=`blastdbcmd -db $la_db -dbtype nucl -entry_batch $ciblestxt`
    run(pipeline(cmd,stdout=extraitsbnk,stderr=devnull))
    vecteurfasta::Vector{String}=[]
    push!(vecteurfasta,fastainit) #seq initiale query

    fileis = open(extraitsbnk) do f
        localfasta::Vector{String}=[]
        localtete::String = ""
        while  !eof(f)
            for l in eachline(f)
                #println(l)
                if startswith(l,'>')
                    #println(l)
                    if ! isempty(localfasta)
                            push!(vecteurtetes,localtete)
                            push!(vecteurtetescoupees,split(strip(localtete,'>'),'=')[1])
                            fantome=join(localfasta,"")
                            #println(localtete,"  ",split(strip(localtete,'>'),'=')[1])
                            push!(vecteurfasta,localtete*'\n'*fantome)
                    end
                    localfasta=[] #initialisation du contenu
                    localtete='>'*split(l,' ')[2]   #affectation de la description
                else 
                    push!(localfasta,strip(l,['\n','-',' ']))
                    
                end
            end
        end
    end


    #println(vecteurtetes)
    collection_avec_query=joinpath(dirutilisateur,"collection_finale.fasta")
    write(collection_avec_query,join(vecteurfasta,'\n'))
    #println(join(vecteurfasta,'\n'))

    return vecteurtetescoupees,vecteurtetes,scoredivers,collection_avec_query

end




function alignement(fastafile)
    return alignefasta(fastafile,"mafft",replace(fastafile,"collection" => "alignement"))
end
function alignefasta(fasta::String,aligneur::String,sorties::String)
    #println("IN   ",fasta," ",aligneur," ",sorties)
    if aligneur == "kalign"
        #println("kalign en route")
        try 
            cmd=`kalign  -n 8 -i $fasta -o $sorties` 
            #println(cmd)
            run(pipeline(cmd,stdout=devnull,stderr=devnull))#,stdout="dev/null",stderr="dev/null"))
            #println(" FAIT    ",fasta," ",aligneur," ",sorties)
        catch
            println(" ERREUR    ",fasta," ",aligneur," ",sorties)
        end
    elseif aligneur == "mafft"
        try
            cmd=`mafft --thread 8 --auto --quiet --amino $fasta`
            run(pipeline(cmd,stdout=sorties,stderr=devnull))#,stdout="dev/null",stderr="dev/null"))
        catch
            println(" ERREUR    ",fasta," ",aligneur," ",sorties)
        end
    elseif aligneur == "famsa"
        try
            cmd=`famsa  -t 1 -v  $fasta  $sorties`
            run(pipeline(cmd,stdout=devnull,stderr=devnull))#,stdout="dev/null",stderr="dev/null"))
        catch
            println(" ERREUR     ",fasta," ",aligneur," ",sorties)
        end
    #cmd = PrgDic["famsa"] + " -t 1 -v " + os.path.join(outDir, "AssQualiteFaite",os.path.basename(FileName).replace("nuc", "Tprot")) + "  " + os.path.join(outDir, "AssQualiteFaite", os.path.basename(FileName).replace("nuc", "prot").replace("fasta", "KTaln"))
    end
    return sorties
end


function trim(a_trimer::String,limitefiltrage::Float64)
    ## pour les tests
    # using SplitApplyCombine 
    # a_trimer="/Users/jean-pierreflandrois/Documents/PKDBGENIESTIPPLE/public/utilisateurs/task_20241007_194634_mXbtzjmP/alignement_finale.fasta"

    # limitefiltrage=0.1
    ##

    listkopf::Vector{String}=[]
    items_arraydespositions::Vector{Vector{Char}}=[]
    matricetransposée::Vector{Vector{Char}}=[]
    
    fileis = open(a_trimer) do f
        localfasta::Vector{String}=[]
        localtete::String = ""
        while  !eof(f) #on utilise la voie rapide de lecture ligne à ligne 
            for l in eachline(f)
                if startswith(l,'>') #on débute un fasta mais...
                    if ! isempty(localfasta) #on vient de lire un fasta entier  donc localfasta est plein et on garde
                        push!(listkopf,localtete) #les noms dans l'ordre de l'alignement
                        push!(items_arraydespositions,collect(join(localfasta,""))) #les matrices position
                    end
                    localfasta=[] #initialisation du contenu
                    localtete=l #split(l,"\n")[2]   #affectation de la description
                else 
                    push!(localfasta,strip(l,['\n']))
                end
            end
        end
    end
    matricetransposée=splitdimsview(combinedimsview(items_arraydespositions,1))
    #on a mantenant les colonnes par position ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'] pour position1 
    # on compte les - 
    matricefiltrage::Vector{Int}=[]
    for u in 1:length(matricetransposée)
        #println(length(matricetransposée))
        if length(filter(x -> x=='-',matricetransposée[u]))/length(matricetransposée[u]) > limitefiltrage #limite des gaps à 10% dans la colonne
            push!(matricefiltrage,u)#le nombre de gaps rapporté a la longueur de la colonne
        end
    end
    #println(sort(matricefiltrage,rev=true)) #on commence par la fin pour garder les positions !!!
    map(x-> deleteat!(matricetransposée,x),sort(matricefiltrage,rev=true))
    # println(length(sort(matricefiltrage,rev=true)))
    # for i in sort(matricefiltrage,rev=true)
    #     println(i)
    #     deleteat!(matricetransposée,i)
    # end
    #on a enlevé les gaps selon la règle
    #####################################################
    # procédure de sélection d'une partie de l'alignement
    # on doit pouvoir réutiliser matricetransposée
    serialize(replace(a_trimer,"alignement_finale.fasta" => "serialisation_inv_mat.ser"), matricetransposée)
    return a_trimer,listkopf,matricetransposée
    #
end

function tranchedemsa(a_trimer::String,debut::Int64,fin::Int64,matricetransposée::Vector{Vector{Char}},listkopf::Vector{String})
    chappe::Bool = true #le télégraphe Chappe :)
    #verifications basiques
    
    
    debut=round((length(matricetransposée)/100)*debut)
    fin=round((length(matricetransposée)/100)*fin)
    dimention=fin-debut
    #println(tranchedemsa,"  ",dimention)
    if dimention < 300
        chappe = false
        return chappe,listkopf,matricetransposée
    end
    #println(length(matricetransposée))
    matricefiltrage::Vector{Int}=[]
    for u in 1:length(matricetransposée)
        #println(u)
        if debut > u 
            push!(matricefiltrage,u)#les colonnes avant min et après max sont délétées
        elseif u > fin 
            push!(matricefiltrage,u)#les colonnes avant min et après max sont délétées
        end
    end
    #println(sort(matricefiltrage,rev=true)) #on commence par la fin pour garder les positions !!!
    #println("finselection")
    map(x-> deleteat!(matricetransposée,x),sort(matricefiltrage,rev=true))
    #println("finmap")
    # println(length(sort(matricefiltrage,rev=true)))
    # for i in sort(matricefiltrage,rev=true)
    #     println(i)
    #     deleteat!(matricetransposée,i)
    # end
    #on a enlevé les gaps selon la règle
    #####################################################
    # procédure de sélection d'une partie de l'alignement
    # on doit pouvoir réutiliser matricetransposée
    serialize(replace(a_trimer,"alignement_finale.fasta" => "serialisation_inv_mat_tranche.ser"), matricetransposée)
    #println("retourpostrimming")
    return chappe,listkopf,matricetransposée #matrice reduite
end



function retourafasta(phase::Int64, a_trimer::String,listkopf::Vector{String},matricetransposée::Vector{Vector{Char}})
    #messages 
    chappe::Bool = true #le télégraphe Chappe :)
    message::String = "ok"
    
    ##################### la matrice transposée est retransposée pour récupérer l'ordre des fasta ######################
    matriceretour = splitdimsview(combinedimsview(matricetransposée,1))
    #on peut vérifier la présence et la position du query et son importance dans l'alignement
    # println(matriceretour[1])
    # println(union(matriceretour[1]))
    query::String=join(matriceretour[1],"") #On prend le query en ligne1
    if phase ==2
        if length(replace(query, '-' => "")) ==0
            chappe = false
            message = "Query is lost, select accurately"
            return chappe, message, ""  #oui tranche == slice :) on ne renvoie rien
        elseif length(replace(query, '-' => "")) < 100
            chappe = false
            message = "Query is really short here, increase it"
            return chappe, message, "" # on ne renvoie rien 
        elseif length(replace(query, '-' => "")) / length(query) < 0.3  
            chappe = false
            message = "Query is not well represented"
            return chappe, message, "" # on ne renvoie rien 
        end
        #subindication=replace("_SelectedRange")
    end
    ################## préparation sorties ######################
    marqueurphase::String=""
    phase == 1 ? marqueurphase = "" : marqueurphase = "_recentrage"
    fastaretour::Vector{String}=[]
    #"# STOCKHOLM 1.0\n"
    stockholmretour::Vector{String}=["# STOCKHOLM 1.0\n"]
    
    posttrimmagefasta=replace(a_trimer,"alignement_finale" => "alignement_assaini"*marqueurphase)
    posttrimmagesth=replace(a_trimer,"alignement_finale.fasta" => "alignement_assaini"*marqueurphase*".sth")

    #on va maintenant enlever les blocs de gaps du query en début et en fin 
    #on pourrait aussi enlever tous les blocs mal alignés
    # debquery=split(query,r"[A-Za-z]")[begin] #blocs --- de debut
    # finquery=split(query,r"[A-Za-z]")[end]     #blocs --- de fin
    # println(length(debquery),"  ",length(finquery))
    # coupure1=length(debquery) >11 ? length(debquery)-10 : length(debquery) #la coupure se fait en amont pour laisser des bouts 
    # coupure2=length(finquery) > 11 ? length(finquery)-10 : length(finquery) #la coupure se fait en aval pour laisser des bouts
    # #on repart dans matricetransposée et on va enlever les blocs n'existant pas dans le query
    # #on commence par la fin pour garder les positions !!!
    # matricefiltrage=sort(union(collect(1:coupure1),collect(coupure2:length(query))),rev=true)
    # map(x-> deleteat!(matricetransposée,x),matricefiltrage) #on enlève
    # matriceretour::Vector{Vector{Char}} = splitdimsview(combinedimsview(matricetransposée,1)) #on remet dans l'ordre des séquences
    
    ############ écriture disque 
    for u in 1:length(matriceretour)
        #println(listkopf[u]*"\n"*join(matriceretour[u],""))
        push!(fastaretour,split(listkopf[u],'=')[1]*"\n"*join(matriceretour[u],""))
        push!(stockholmretour,replace(strip(split(listkopf[u],'=')[1],'>'),' '=> "")*"    "*join(matriceretour[u],""))
    end
    finalementfasta::String=join(fastaretour,"\n")
    write(posttrimmagefasta,finalementfasta)
    finalementstockholm::String=join(stockholmretour,"\n")*"\n//"
    write(posttrimmagesth,finalementstockholm)
    return (chappe, message, posttrimmagesth)
end

function lis_moi_fasta(entree::String)
    listkopf::Vector{String}=[]
    vecteurdessequences::Vector{String}=[]
    fileis = open(entree) do f
        localfasta::Vector{String}=[]
        localtete::String = ""
        while  !eof(f) #on utilise la voie rapide de lecture ligne à ligne 
            for l in eachline(f)
                if startswith(l,'>') #on débute un fasta mais...
                    if ! isempty(localfasta) #on vient de lire un fasta entier  donc localfasta est plein et on garde
                        push!(listkopf,localtete) #les noms dans l'ordre de l'alignement
                        push!(vecteurdessequences,join(localfasta,""))
                    end
                    localfasta=[] #initialisation du contenu
                    localtete=l #split(l,"\n")[2]   #affectation de la description
                else 
                    push!(localfasta,strip(l,['\n']))
                end
            end
        end
    end
    return vecteurdessequences
end
function panoramatographe_nuc(entree::String,sortie::String,cotécarré::Int) #png ici
    A::Vector{String}=lis_moi_fasta(entree)
    if occursin("collection_finale.fasta",entree) #ajout des "-" par rapport à la sequence la plus longue
        lmax=sort(map(x -> length(x),A),rev=true)[1]
        A=map(x-> x*"-"^(lmax-length(x)),A)
    end
    # a partir de A on fabrique une liste de colorisation et une liste de nom de bases
    B::Vector{Vector{SubString{String}}}= [split(i,"") for i in A] 
    colorisé::Vector{Vector{String}}=map((j) -> [replace(i,"A" => "red", "T" => "blue","C" => "green","G" => "yellow","N" => "grey","-" => "black") for i in j],B)
    couleursvector::Vector{String}= [(colorisé...)...] #la liste de 1 ... n
    bases::Vector{Char}=[(A...)...] # les noms des bases de 1 ... n
    type::String = "png"
    sortie=sortie*"."*type
    #println(length(B[1])*cotécarré,"  " ,length(A)*cotécarré)
    @png begin #defaut  pour le néobibi avec un carré de 15 ou 20 pour lettres sinon 5
        background("black")
        t::Table = Table(length(B), length(B[1]), cotécarré,cotécarré) # rows, columns, wide, high
        for (pt, n) in t
            setcolor(couleursvector[n])
            box(pt, cotécarré, cotécarré, :fill)
            setcolor("black")
            # if cotécarré >=10
            #     fontsize(cotécarré/2)
            #     text(string(bases[n]), t[n], halign=:center, valign=:middle)
            # end
        end
    end length(B[1])*cotécarré length(A)*cotécarré sortie #la taille de l'alignement
    return sortie
end

# function faitrapidnj(pourarbre::String)
#     sorties=replace(pourarbre,"sth" => "ph")
#     cmd=`rapidnj $pourarbre  -i sth -t d -a kim -o t -c 8 -b 100 -n -x $sorties`
#     run(pipeline(cmd,stdout=sorties,stderr=devnull))
#     raciné=replace(pourarbre,"sth" => "midpoint.ph")
#     # if isfile("/usr/local/bin/fastroot.py")
#     #     cmd=`python3.11 /usr/local/bin/fastroot.py -i $sorties -o $raciné -m MP`
#     # else
#     #     cmd=`python3.11 /Library/Frameworks/Python.framework/Versions/3.11/bin/fastroot.py -i $sorties -o $raciné -m MP` #cas de la maison
#     # end
#     run(pipeline(cmd,stdout=raciné,stderr=devnull))
#     baum=replace(raciné,".midpoint.ph" => ".tree.pdf")
#     cmd=`Rscript public/R/baumphangorn.R $sorties $raciné $baum` 
#     #println(cmd)
#     run(pipeline(cmd,stdout=devnull,stderr=devnull))#,stdout="dev/null",stderr="dev/null"))
#     return baum
# end
function faitrapidnj(pourarbre::String) #format stockholm 
    sorties=replace(pourarbre,"sth" => "ph")
    cmd=`rapidnj $pourarbre  -i sth -t d -a kim -o t -c 8 -b 100 -n -x $sorties`
    #println(cmd)
    run(pipeline(cmd,stdout=sorties,stderr=devnull))
    #nouveaux noms
    raciné=replace(pourarbre,"sth" => "midpoint.ph")
    baum=replace(raciné,".midpoint.ph" => ".tree.svg") #
    #envoi R
    cmd=`Rscript R/baumphangorn.r $sorties $raciné $baum` 
    #println(cmd)
    run(pipeline(cmd,stdout=devnull,stderr=devnull))
    return baum
end


function exemples(demande::String) #réaction à une demande de démonstration 
    if demande == "trimming1" 
        return ">trimming_exemple\nGGGAATCTTAGACAATGGGGGCAACCCTGATCTAGCCATGCCGCGTGAGTGATGAAGGCCCTAGGGTTGTAAAGCTCTTTCAGCTGGGAAGATAATGACGGTACCAGCAGAAGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGACCTGAAAGTTGGGGGTGAAATCCCGGGGCTCAACCTCGGAACTGCCTTCAAAACTACGGGTCTGGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATGCCAGTCGTCGGGCAGCATGCTGTTCGGTGACACACCTAACGGATTAAGCATTCCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAACCCTTGACATCGCAGGACCGCCCGAGAGATCGGGTTTTCTCGTAAGAGACCTGTGGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTCGGTTAAGTCCGGCAACGAGCGCAACCCACGTCCCTAGTTGCCAGCATTCAGTTGGGCACTCTATGGAAACTGCCGATGATAAGTCGGAGGAAGGTGTGGATGACGTCAAGTCCTCATGGCCCTTACGGGTTGGGCTACACACGTGCTACAATGGTGGTGACAGTGGGTTAATCCCCAAAAGCCATCTCAGTTCGGATTGGGGTCTGCAACTCGACCCCATGAAGTTGGAATCGCTAGTAATCGCGGAACAGC"
    else 
    return("pas trouvé")
    end
end

"""
Exemple MINIMAL qui vérifie toutes les fonctions :)

SEQ=">gi|646115914|gb|KJ479307.1| Uncultured bacterium clone A5521f07 16S ribosomal RNA gene, partial sequence\nTGGGGAGCAAACAGGATTAGATACCCTGGTAGTACCACGCTGTAAACGATGGGTACTAGGTGTAGGGGGC\nGTAGCCCTCTGTGCCGCAGTTAACACATTAAGTACCCCGCGTGGGAAGTACGATCGCAAGATTAAAACTC\nAAAGGAATTGACGGGGGCCCGCACAAGCAGCGGAGCATGTGGTTTAATTAGAAGCAACGCGAAGAACCTT\nACCAGGGCTTGACATCGTACTAACCCTGTGGAAACACGGGGGTGCCCTTCGGGGAACGGTGAGACAGGTG\nGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCACAACGAGCGCAACCCTTGCCT\nTTAGTTGCCAGCAATTCGGTTGGGCACTCTAGAGGGACTGCCTGGGTTAACCAGGAGGAAGGTGGGGATG\nACGTCAAATCATCATGCCCCTTATGCTCTGGGCTACACACGTGCTACAATGGCCGGTACAATGAGTTGCA\nAACCCGCGAGGGGGAGCTAATCTCAAAAACCGGTCCCAGTTCGGATTGTAGGCTGCAACTCGCCTACATG\nAAGCTGGAGTTGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCG\nCCCGTCACACCATGAGAGCCGGCAACACCCGAAGCCAGTGGGCTAACCCGCAAGGGAGGCAGCTGTCGAA\nGGTGGGGTTGGTGATTG"
p = replace(SEQ,'>' => ">QUERY_"  ," " => "_", "," => "_", ";" => "_", ":" => "_", "(" => "_", ")" => "_", "__" => "_")
vecteurtetescoupees,vecteurtetes,scoredivers,collection_avec_query = faitblast(bnk,p,diris,30)


alignementfait = alignement("/Users/jean-pierreflandrois/Documents/PkPhyExplo/public/utilisateurs/task_20241019_124721_rM2TmZbG/atelier_20241019_124721_rM2TmZbG/collection_finale.fasta")

a_trimer,listkopf,Mat = trim("/Users/jean-pierreflandrois/Documents/PkPhyExplo/public/utilisateurs/task_20241019_124721_rM2TmZbG/atelier_20241019_124721_rM2TmZbG/alignement_finale.fasta",0.1)
matricetransposée = copy(Mat)
chappe, message, posttrimmage = retourafasta(1,"/Users/jean-pierreflandrois/Documents/PkPhyExplo/public/utilisateurs/task_20241019_124721_rM2TmZbG/atelier_20241019_124721_rM2TmZbG/alignement_finale.fasta",listkopf,matricetransposée)
panoramatographe_nuc(replace(posttrimmage,"sth" => "fasta"),replace(posttrimmage,"sth" => "image"),4)
faitrapidnj(posttrimmage)
matricetransposée = copy(Mat)

chappe2,listkopf2,matricetransposée2 = tranchedemsa(a_trimer, 800, 1600, matricetransposée, listkopf)
chappe, message, posttrimmage = retourafasta(2,"/Users/jean-pierreflandrois/Documents/PkPhyExplo/public/utilisateurs/task_20241019_124721_rM2TmZbG/atelier_20241019_124721_rM2TmZbG/alignement_finale.fasta",listkopf2,matricetransposée2)
panoramatographe_nuc(replace(posttrimmage,"sth" => "fasta"),replace(posttrimmage,"sth" => "image"),4)
faitrapidnj(posttrimmage)

matricetransposée = copy(Mat)
println(length(matricetransposée))
chappe3,listkopf3,matricetransposée3 = tranchedemsa(a_trimer, 800, 1400, matricetransposée, listkopf)
chappe4, message4, posttrimmage4 = retourafasta(2,"/Users/jean-pierreflandrois/Documents/PkPhyExplo/public/utilisateurs/task_20241019_124721_rM2TmZbG/atelier_20241019_124721_rM2TmZbG/alignement_finale.fasta",listkopf3,matricetransposée3)
panoramatographe_nuc(replace(posttrimmage4,"sth" => "fasta"),replace(posttrimmage4,"sth" => "image"),4)
faitrapidnj(posttrimmage)
"""


end #Module end
