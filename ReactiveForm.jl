module ReactiveForm
using GenieFramework
using Stipple
using StippleDownloads
using .Main.App.Analysis
using DataFrames
@genietools

@app begin
    @in trigger = false
    @in clearit = false
    @in choixposttrim = false
    @in download_event = false
    @in montre_moi_tirer = false
    @in blast_fait= false
    @in collect_fait= false
    @in align_fait= false
    @in trim_fait= false
    @in trim_fait_persistant= false
    @in arbre_fait= false
    @in figure_arbre_fait= false
    @in selectionfaite =false
    @in banqueblast = "TRECS_16SrRNA.fst"
    @in termine = "Blast running"
    @in banqueselectionnée = "TRECS_16SrRNA.fst"
    # @out listofbanques = ""
    @in listebanques = ["TRECS_16SrRNA.fst","TRECS_23SrRNA.fst","TRECS_5SrRNA.fst"]
    @in borne_longueur_inf=800
    @in borne_longueur_sup=2000
    @in requestedseq = 50
    @in limitedsearch = 50
    @in selectionintervalle = RangeData(0:100)
    @out selectionmsainf=0
    @out selectionmsasup=100
    @out drapeau_tranche::Bool = false
    ############################
    # 0:5:100,
    # :selectionintervalle,
    # markers = true,
    # var:"marker-labels" = true,
    # color = "secondary"

    # @in selectionintervalle::RangeData{Int}=RangeData(0:100)
    # @out Range_Markers_labels =
    #     [Dict(:value => i, :label => string(i) * "%") for i = 0:10:100]
    # @onchange selectionintervalle
    #     println(selectionintervalle, typeof(selectionintervalle))
    # end
        #@in paramètres = false
    ############################
    @onchange banqueselectionnée begin
        banqueblast = banqueselectionnée
        if banqueblast == "TRECS_16SrRNA.fst"
            borne_longueur_inf=800
            borne_longueur_sup=2200
        elseif banqueblast == "TRECS_23SrRNA.fst"
            borne_longueur_inf=1400
            borne_longueur_sup=4500
        elseif banqueblast == "TRECS_5SrRNA.fst"
            borne_longueur_inf=75
            borne_longueur_sup=200
        end
        
        println("selection :",borne_longueur_inf," ",borne_longueur_sup)
        
    end
    @onchange requestedseq begin
        limitedsearch = requestedseq
        #println(requestedseq," help seq",  banqueblast,"  ",limitedsearch)
        # println("RM :",requestedseq," ",typeof(requestedseq))
        #println("selection :",limitedsearch," ",typeof(limitedsearch))
        #paramètres = true
    end
    
    @out travail = false
    @out travail2 = false
    @in ddff_pagination = DataTablePagination(rows_per_page = 1)

    @out m = "non actif "
    @in S = ""
    @out p = ""
    @out resu= ()
    #@out BL=([],[],[],"")

    @in posdsk::String ="waiting for the link"
    @in downloadinfo::String = "not Ready"
    #@out genomeid = []
    @out genomeid::Vector{String}=[]
    @out vecteurtetes::Vector{String}=[]
    @out scoredivers::Vector{Tuple{SubString{String}, SubString{String}}}=[] # genomeid,vecteurtetes,scoredivers,fichierfasta
    ##ajoutscar :  vecteurtetescoupees,lataxinomie,evalue,scores,ali_length,identitynumber,identitypc,gapsopen,collection_avec_query
    @out lataxinomie::Vector{String}=[]
    @out qualité::Vector{String}=[]
    @out ali_length::Vector{Int64}=[]
    @out identitynumber::Vector{Int64}=[]
    @out identitypc::Vector{Float64}=[]
    @out gapsopen::Vector{Int64}=[]

    @out fichierfasta::String = ""
    @out sscores::Vector{Int64}=[]
    @out evalue::Vector{Float64}=[]

    #@out fichierfasta = ""
    @out estalign = ""
    @out fasta_trimé =""
    @out fintrim = ""
    @out chappe::Bool = false
    @out message = ""
    @out termine = "Ready for a new submission"
    #listedesaextraire,lataxinomie,laqualite,evalue,scores,ali_length,identitynumber,identitypc,gapsopen,collection_avec_query
    @out ddff = DataTable(DataFrame(SpeciesGenome=String["no data101"],GenomeId=String["no data"],Quality=String["no data"],Scores=String["NAN"],Evalue=String["NAN"],Ali_Length=String["NAN"],IdentityNumber=String["NAN"],Identity=String["NAN"],GapsOpen=String["NAN"]))# sseqid           stitle                             evalue   bitscore  length  nident  pident   gapopen
    @out matricetransposee::Vector{Vector{Char}} =[] #la matrice transposée initiale après trimming
    @out transposée_msa::Vector{Vector{Char}} =[]
    @out matricetranchetransposee::Vector{Vector{Char}} =[] #la matrice transposée en cas de selection d'une tranche dans le MSA
    @out listkopf::Vector{String}=[]
    @out baumist = ""
    @out baumist_sel = ""
    @out seaview_a = ""
    @out seaview_b = ""
    @out seaview_c = ""
    @out seaview_sel = ""
    @out borne_longueur_inf = 600 #dépend de la banque ici 16S
    @out borne_longueur_sup = 2500
    @out pourgzip = ""
    
    @onbutton clearit begin
        trigger = false
        iclearit = false
        download_event = false
        montre_moi_tirer = false
        blast_fait= false
        collect_fait= false
        align_fait= false
        trim_fait= false
        trim_fait_persistant= false
        arbre_fait= false
        figure_arbre_fait= false
        selectionfaite =false
        travail = false
        travail2 = false
        ddff_pagination = DataTablePagination(rows_per_page = 1)

        m = "non actif "
        S = ""
        p = ""
        resu= ()
        BL=([],[],[],"")
        posdsk ="waiting for the link"
        downloadinfo = "not Ready"
        genomeid = []
        nomslongs = []
        scores = []
        evalues = []
        fichierfasta = ""
        
        estalign = ""
        fintrim = ""
        termine = "Cleared, Ready for a new submission, Verify the Blast Bank and the number of selected hits"
        ddff = DataTable(DataFrame(SpeciesGenome=String["no data157"],GenomeId=String["no data"],Quality=String["no data"],Scores=String["NAN"],Evalue=String["NAN"],Ali_Length=String["NAN"],IdentityNumber=String["NAN"],Identity=String["NAN"],GapsOpen=String["NAN"]))
        baumist = ""
        seaview_a = ""
        seaview_b = ""
        seaview_c = ""
        seaview_sel = ""
        baumist_sel = ""
        # borne_longueur_inf = 600
        # borne_longueur_sup = 3500
        pourgzip = ""
        matricetrimtransposée::Vector{Vector{Char}}=[]
        
    end
    @onbutton trigger begin
        #nettoyer avant en cas de non utilisation du bouton clear
        ddff = DataTable(DataFrame(SpeciesGenome=String["no data172"],GenomeId=String["no data"],Quality=String["no data"],Scores=String["NAN"],Evalue=String["NAN"],Ali_Length=String["NAN"],IdentityNumber=String["NAN"],Identity=String["NAN"],GapsOpen=String["NAN"]))
        posdsk ="waiting for the link"
        downloadinfo = "not Ready"
        genomeid = []
        nomslongs = []
        scores = []
        evalues = []
        fichierfasta = ""
        estalign = ""
        fintrim = ""
        baumist = ""
        baumist_sel = ""
        seaview_a = ""
        seaview_b = ""
        seaview_c = ""
        seaview_sel = ""
        pourgzip = ""
        download_event = false
        montre_moi_tirer = false
        blast_fait= false
        collect_fait= false
        align_fait= false
        trim_fait= false
        trim_fait_persistant= false
        selectionfaite =false
        arbre_fait= false
        figure_arbre_fait= false
        #fin nettoyage
        S=strip(S)
        ddff_pagination = DataTablePagination(rows_per_page = 5)
        # if paramètres == false
        #     termine = "Please choose and valid the parameters"
        termine = "Fasta verifications"
        if S == ""                      #pas vide
            termine = "MISSING FASTA, PLEASE SUBMIT A SEQUENCE IN FASTA FORMAT\nsubmit ? here for help"
            S = "MISSING FASTA ! "
        elseif occursin("?length",S)
            termine = "See help in the fasta box"
            S = "A decent phylogenetic reconstruction needs suffisant phylogenetic signal\n and a balanced phylogenetic signal with that of the sequences in the database\n With $banqueblast min=$borne_longueur_inf, max=$borne_longueur_sup bp"

        elseif occursin("?demo:",S)
            termine = "Selection of a demonstrative sequence"
            println(strip(split(S,"demo:")[2]))
            S = exemples(strip(split(S,"demo:")[2]))
            println(S)
            if S == "pas trouvé"
                termine = "No corresponding demonstrative sequence"
                S = "?demo:trimming1 "#...
            end
        elseif occursin('?',S)
            termine = "See help in the fasta box"
            S = "The PK Phylogeny Explorer works on nucleic sequences only\nFasta format is :\n>sequence_information\nATCGGCT....TNTCGGATT\n\n"
        elseif validate_fasta_nucleic(S,true) == false
            if validate_fasta_proteic(S,true) == true
                termine = "SEEMS TO BE A PROTEIN SEQUENCE, ONLY NUCLEIC HERE"
                S = " SEEMS TO BE A PROTEIN SEQUENCE, ONLY NUCLEIC ALLOWED\nsubmit ? here for help"
            else
                termine = "ILL FORMATED FASTA OR NOT CLEARED, NOTE THAT MULTI-FASTA SUBMISSION IS NOT ALLOWED\nsubmit ? here for help"
                S = " ILL FORMATED FASTA OR FORM NOT CLEARED "
            end
        
        else #@private a essayer
            
            longueurdefasta=length(replace(split(S,"\n",limit=2)[2],"\n" => "","-" => ""))
            decisonlongueur_l= longueurdefasta < borne_longueur_inf ? false : true 
            decisonlongueur_L= longueurdefasta > borne_longueur_sup ? false : true 
            if ! decisonlongueur_L
                termine = "A CANDIDATE SEQUENCE FOR A BLAST SEARCH AGAINST THE DB MUST BE SHORTER (max bp $borne_longueur_sup)"
                S = "help: enter  ?length in the fasta form"
            else
                fairearbre= true
                if decisonlongueur_l == false
                    termine = "TREE WILL NOT BE BUILT: BLAST ONLY  (min bp $borne_longueur_inf)"
                    #S = " THE LENGTH OF THE QUERY IS INADEQUATE FOR PHYLOGENY BLAST ONLY"
                    #temporisation ici 
                    sleep(2)
                    fairearbre= false
                end
                #pour tester et ne pas bloquer :
                #S=">PROBLEM !! THIS IS A DEMO\nCACGCTGGAATGCCGGGACCCACAAACGCTCCGGCGCTGCAGGATGCGGCTGCGGCCGATTAGGTAGACGGTGGGGTAACGGCCCACCGTGCCGATAATCGGTACGGGTTGTGAGAGCAAGAGCCCGGAGACGGAATCTGAGACAAGATTCCGGGCCCTACGGGGCGCAGCAGGCGCGAAACCTTTACACTGCACGCAAGTGCGATAAGGGGACTCCGAGTGCGAGGGCATATAGTCCTCGCTTTTGTGAATCGTAAGGTGATTCACGAATAAGAGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAATACCGGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCCGGGGTAGGAGTGAAATCCTGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAACCGGATTAGATACCCGGGTAGTCCTAGCCGTAAACAATGTTCACTAGGTGTGGCACAGGCTACGAGCCTGTGCTGTGCCGTAGGGAAGCCGAGAAGTGAACCGCCTGGGAAGTACGTCCGCAAGGATGAAACTTAAAGGAATTGGCGGGGGAGCACTACAACCGGAGGAGCCTGCGGTTTAATTGGACTCAACGCCGGACATCTCACCAGCTCCGACTACAGTGATGACGATCAGGTTGATGACCTTATCACGACGCTGTAGAGAGGAG"
                p = replace(S,'>' => ">QUERY_"  ," " => "_", "," => "_", ";" => "_", ":" => "_", "(" => "_", ")" => "_", "__" => "_")
                m = "pressed"
                travail = true
                termine = "Creation of the directory"
                #p = replace(S,'>'=>">QUERY_")
                
                # m = length(N)
                
                #println(p)
                
                posdsk = uniqueutilisateur()
                termine = "Blast running and collecting similar sequences"
                #println(typeof(p),p," ** ** ", typeof(posdsk),posdsk)
                genomeid,lataxinomie,qualité,evalues,scores,ali_length,identitynumber,identitypc,gapsopen,fichierfasta = faitblast(banqueblast,p,posdsk,limitedsearch) ##listedesaextraire,lataxinomie,laqualite,evalue,scores,ali_length,identitynumber,identitypc,gapsopen,collection_avec_query
                # BL = faitblast(banqueblast,p,posdsk,limitedsearch)
                # if BL[1] ≠ [] #No hits found
                #     fichierfasta=BL[4]
                #     genomeid = BL[1]

                #     for u in 1:length(genomeid)
                #         push!(scores,BL[3][u][1])
                #         push!(evalues,BL[3][u][2])
                #     end
                if genomeid ≠ [] #No hits found
                    # fichierfasta=collection_avec_query
                    # genomeid = vecteurtetescoupees
                    #rappel SpeciesGenome=String["no data"],Scores=String["NAN"],Evalue=String["NAN"],Ali_Length=String["NAN"],Nb_Identity=String["NAN"],Identity=String["NAN"],GapsOpen=String["NAN"]
                    ddff = DataTable(DataFrame(SpeciesGenome=lataxinomie,GenomeId=genomeid,Quality=qualité,Scores=scores,Evalue=evalues,Ali_length=ali_length,Nb_Identity=identitynumber,Identity=identitypc,GapsOpen=gapsopen))# rempplir DataTable(DataFrame(SpeciesGenome=String["no data"],Scores=String["NAN"],Evalue=String["NAN"],ali_length=String["NAN"],identitynumber=String["NAN"],identitypc=String["NAN"],gapsopen=String["NAN"]))
                    blast_fait = true
                  
                    if fairearbre #on fait l'alignement et l'arbre initial avec le trimming 0.1 habituel
                        seaview_a = panoramatographe_nuc(fichierfasta,replace(fichierfasta,".fasta" =>""),1)#(entree::String,sortie::String,cotécarré::Int)
                        seaview_a = split(seaview_a,"public/")[2]
                        collect_fait = true
                        termine = "Alignment running"
                        estalign = alignement(fichierfasta)
                        align_fait = true
                        termine = "Alignment viewer"
                        seaview_b = panoramatographe_nuc(estalign,replace(estalign,".fasta" =>""),1)#(entree::String,sortie::String,cotécarré::Int)
                        seaview_b = split(seaview_b,"public/")[2]
                        #println(seaview_b)
                        termine = "Trimming running"
                        #println("306",estalign)
                        fasta_trimé,listkopf,transposée_msa = trim(estalign,0.1) #la fonction est coupée en deux pour permettre le greffage de la sélection de tranche
                        #println("308",fasta_trimé)

                        matricetrimtransposée = copy(transposée_msa) #important pour éviter les modifications NB matrice après trimming seulement (on détruit la matrice initiale)
                        chappe, message, fintrim = retourafasta(1,fasta_trimé,listkopf,matricetrimtransposée) #ancienne fonction etait directe mais il nous faut la matrice transposée pour les tranches à venir
                        #println(chappe,"***",message)
                        trim_fait = true
                        trim_fait_persistant= true
                        seaview_c = panoramatographe_nuc(replace(fintrim,".sth" =>".fasta"),replace(fintrim,".sth" =>""),1)
                        seaview_c = split(seaview_c,"public/")[2]

                        #p("retour trim $fintrim")
                        termine = "Tree constuction"
                        rapidnj=faitrapidnj(fintrim)
                        baumist=split(rapidnj,"public/")[2]
                        figure_arbre_fait = true
                        termine = " Tree available "#*baumist
                        montre_moi_tirer = true
                        # pourgzip=compresser(splitdir(posdsk)[1]) #compresser atelier sous le nom de l'utilisateur
                        # println("POUR GZIP",pourgzip)
                    else
                        termine = "THE LENGTH OF THE QUERY IS UNADEQUATE FOR PHYLOGENY: BLAST ONLY. For help: enter  ?length in the fasta form"
                    end
                else
                    termine = "***** No hits found ***** Did you choose the right bank ?"

                end
                #reininitialisation
                S = "EXPLORATION FINISHED\n-------------------\n"*S*"\n-------------------\nEXPLORATION FINISHED\n-------------------\n"
                #println("*************")
                #enlevé soir
                travail = false
                # #println(S)
                trigger = false
                clearit = false
                # m = "non actif "
                # p = ""
                # resu= ()
                # BL=([],[],[],"")
                # #posdsk =""
                # genomeid = []
                # nomslongs = []
                # scores = []
                # evalues = []
                # fichierfasta = ""
                # estalign = ""
                # fintrim = ""
                # #seaview = ""
                # limitedsearch= 50
                #ddff = DataTable(DataFrame(SpeciesGenome=String["no data"],Scores=String["NAN"],Evalue=String["NAN"]))
                #baumist = ""
                #println("*************")
                #uploader(label="Upload Image", autoupload=true, multiple=true, method="POST", url="/upload", field__name="img")
                #/Users/jean-pierreflandrois/Documents/PKDBGENIESTIPPLE/public/utilisateurs/task_20240929_123919_kLVW5Yde
                #Router.serve_static_file(joinpath(split(posdsk,"public/")[2],"blast.out.html")) #public/Users/jean-pierreflandrois/Documents/PKDBGENIESTIPPLE/public/utilisateurs/task_20240929_133908_PVboKQVF/alignement_assaini.tree.pdf
                # Router.serve_static_file("public",joinpath(split(posdsk,"public/")[2],"alignement_assaini.tree.pdf"))
            end
        end
    end
    @event download_event begin
        downloadinfo = "running..." #/Users/jean-pierreflandrois/Documents/PKDBGENIESTIPPLE/public/utilisateurs/task_20241013_220057_lfXiwPpZ/atelier
        #println("downloda process")
        #println(posdsk)
        pourgzip=compresser(splitdir(posdsk)[1]) #compresser atelier sous le nom de l'utilisateur
        #println("=====")
        #println(pourgzip)
        downloadinfo = pourgzip
        try
            solution = joinpath("public", "utilisateurs", splitdir(posdsk)[1],pourgzip)
            #println("downloading...",joinpath("public", "utilisateurs", splitdir(posdsk)[1],pourgzip))
            io = IOBuffer()
            open(solution, "r") do file 
                write(io, read(file))
            end
            seekstart(io)
            download_binary(__model__,take!(io), pourgzip)
            downloadinfo = "Done !"
            #on va maintenant vider l'écran
            trigger = false
            iclearit = false
            download_event = false
            montre_moi_tirer = false
            blast_fait= false
            collect_fait= false
            align_fait= false
            trim_fait= false
            trim_fait_persistant= false
            selectionfaite =false
            arbre_fait= false
            figure_arbre_fait= false
                #banqueblast = ""
                #banqueselectionnée="TRECS_16SrRNA.fst"
                    # @out listofbanques = ""
                #listebanques = ["TRECS_16SrRNA.fst","TRECS_23SrRNA.fst"]
                #choixbanque = false
                #requestedseq = 50
                #limitedsearch = 50
                #paramètres = false
            travail = false
        
            ddff_pagination = DataTablePagination(rows_per_page = 1)

            m = "non actif "
            S = ""
            p = ""
            resu= ()
            BL=([],[],[],"")
            
            downloadinfo = "not Ready"
            genomeid = []
            nomslongs = []
            scores = []
            evalues = []
            fichierfasta = ""
            
            estalign = ""
            fintrim = ""
            termine = "Cleared, Ready for a new submission, Verify the Blast Bank and the number of selected hits"
            ddff = DataTable(DataFrame(SpeciesGenome=lataxinomie,GenomeId=genomeid,Quality=qualité,Scores=scores,Evalue=evalues,Ali_Length=ali_length,Nb_identity=IdentityNumber,Identity=identitypc,GapsOpen=gapsopen))
            #Ali_Length=String["NAN"],IdentityNumber=String["NAN"],Identity=String["NAN"],GapsOpen=String["NAN"]
             #genomeid,lataxinomie,qualité,evalues,scores,ali_length,identitynumber,identitypc,gapsopen,fichierfasta
            baumist = ""
            seaview_a = ""
            seaview_b = ""
            seaview_c = ""
            borne_longueur_inf = 600
            borne_longueur_sup = 3500
            pourgzip = ""
            matricetrimtransposée::Vector{Vector{Char}}=[]
            ################  et maintenant détruire le classeur utilisateur #################
            #rm(posdsk, recursive=true, force=true)
            #vider l'information correspondante
            posdsk ="waiting for the link"
        catch ex
            println("Error during download: ",ex)
        end
    end

    @event choixposttrim begin
        #println("dans la selection posttrim")
        #travail = true
        travail2 = true
        selectionmsainf::Int64 = selectionintervalle.range.start
        selectionmsasup::Int64 = selectionintervalle.range.stop
        drapeau_tranche = true
        #println("cool de ",selectionmsainf," a ",selectionmsasup,"  ",typeof(selectionmsainf))
        trim_fait = false
        trim_fait_persistant= true
        
        matricetrimtransposée = copy(transposée_msa) #important pour éviter les modifications NB matrice après trimming seulement (on détruit la matrice initiale)
        #println("copiefaite",selectionmsainf,"   ",selectionmsasup,"  ",listkopf)
        #println(matricetrimtransposée)
        #chappe2,listkopf2,matricetransposée2 = tranchedemsa(a_trimer, 800, 1600, matricetransposée, listkopf)
        chappe,listkopf,transposée_msa2 = tranchedemsa(fasta_trimé,selectionmsainf,selectionmsasup,matricetrimtransposée,listkopf)#listkopf,matricetrimtransposée
        #println("fait 400","  ",chappe)
        #chappe, message, posttrimmage = retourafasta(2,"/Users/jean-pierreflandrois/Documents/PkPhyExplo/public/utilisateurs/task_20241019_124721_rM2TmZbG/atelier_20241019_124721_rM2TmZbG/alignement_finale.fasta",listkopf2,matricetransposée2)
        chappe, message, fintrim = retourafasta(2,fasta_trimé,listkopf,transposée_msa2)#ancienne fonction etait directe mais il nous faut la matrice transposée pour les tranches à venir
        #println("fait 403"," ",chappe,"  ",fintrim)
        #println(chappe,"***",message)
        # trim_fait = true
        selectionfaite =true
        #println("--- FIN ---")
        seaview_sel = panoramatographe_nuc(replace(fintrim,".sth" =>".fasta"),replace(fintrim,".sth" =>""),1)
        seaview_sel = split(seaview_sel,"public/")[2]
        #chappe2,listkopf2,matricetransposée2 = tranchedemsa(a_trimer, 800, 1600, matricetransposée, listkopf)
        #chappe, message, posttrimmage = retourafasta(2,"/Users/jean-pierreflandrois/Documents/PkPhyExplo/public/utilisateurs/task_20241019_124721_rM2TmZbG/atelier_20241019_124721_rM2TmZbG/alignement_finale.fasta",listkopf2,matricetransposée2)
        #panoramatographe_nuc(replace(posttrimmage,"sth" => "fasta"),replace(posttrimmage,"sth" => "image"),4)
        #faitrapidnj(posttrimmage)
        
        #p("retour trim $fintrim")
        termine = "Tree constuction"
        rapidnj=faitrapidnj(fintrim)
        baumist_sel=split(rapidnj,"public/")[2]
        figure_arbre_fait = true
        termine = " Tree available "#*baumist_sel
        montre_moi_tirer = true
        #travail = false
        travail2 = false
    end 

end
#[btn("Trigger action", @click(:trigger))]#
function ui()  #btn("valider",color="red",@click("press_btn = true")), # @onbutton start begin
    # cell([
    #     textfield(type="textarea","fasta ?", :N),btn("valider",color="red",@click("press_btn = true")),
    #     p("fasta est {{N}} de longueur is {{m}}") # random numbers is {{m}}
    # ])
    #[btn("Trigger action", @click(:trigger)),p("{{m}}")]
    [cell([h1("Prokaryotes Phylogenetics Explorer [nuc] Workshop")])#toolbar("Configuration", class = "bg-primary text-white shadow-2")
    cell([h6("LBBE UMR5558 Université Lyon1-CNRS & Master bioinfo@lyon Université Lyon1"),  a(href="https://github.com/jpflandrs/PkXplore","Code")])
    separator(color = "primary")
    p(h5("Global Parameters"))
    #cell([h5("Global Parameters")])
    row([
        column(class="q-pa-sm",Stipple.select(:banqueselectionnée, options=:listebanques, emitvalue=true, clearable=true, useinput=true, counter = true, fillinput=true, filled = true, label="Bank Selection"), sm=6),
        column(class="q-pa-sm",cell([p("Nb seqs"),slider(30:5:75, :requestedseq, var"marker-labels" = true, color = "grey")]), sm=6),
        #column(class="q-pa-sm",btn("Validate", @click(:choixbanque)), sm=1),
        ])
    cell([separator(color = "primary")])
        
    cell(h5("Query"),[textfield(type="textarea","sequence fasta ?", :S),btn("Submit", @click(:trigger),loading =:travail,color = "secondary"),btn("Clear", @click(:clearit))])
        
    cell(["Process info: {{termine}}"])
        Html.div(cell([h5("Blast against DB:{{banqueselectionnée}}"),table(:ddff, paginationsync = :ddff_pagination, flat = true, bordered = true, title = "BLAST")]), @showif("blast_fait") )
        Html.div(cell([h5("Collected sequences"), imageview(
            src = :seaview_a,
            style = "width: 100%"
        )]),@showif("collect_fait"))
    Html.div(cell([h5("Aligned sequences"), imageview(
            src = :seaview_b,
            style = "width: 100%"
        )]), @showif("align_fait"))
    Html.div(cell([h5("Trimmed Alignment"), imageview(
            src = :seaview_c,
            style = "width: 100%"
        )]), @showif("trim_fait_persistant"))
    
        ######################################
    Html.div(cell([h5("Select a zone and redo"), Stipple.range(
            0:1:100,
            :selectionintervalle,
            markers = true,
            #var:"marker-labels" = true,
            color = raw"secondary"
            ) ,cell(["Tranche info: {{selectionmsainf}}% - {{selectionmsasup}}%"]), 
            #btn("Select Slice", @click(:choixposttrim), loading = :travail2, color = "secondary")]), 
            btn("Select", @on(:click,:choixposttrim), loading = :travail2, color = "secondary")]),
            @showif("trim_fait"))
        #########################################
    Html.div(cell([h5("Selected Slice of Trimmed Alignment"), h6("Slice {{selectionmsainf}}% - {{selectionmsasup}}%"), imageview(
            src = :seaview_sel,
            style = "width: 100%"
        )]), @showif("selectionfaite"))
    Html.div(cell([h5("Whole Trimmed Alignment : Phylogenetic reconstruction"), imageview(
            src = :baumist,
            style = "width: 100%" 
        )]),@showif("figure_arbre_fait"))
    Html.div(cell([h5("Selected Slice of Trimmed Alignment : Phylogenetic reconstruction"), h6("Slice {{selectionmsainf}}% - {{selectionmsasup}}%"), imageview(
            src = :baumist_sel,
            style = "width: 100%"
        )]),@showif("selectionfaite"))
        
        cell([separator(color = "primary")])
    Html.div(
        row([cell([ btn(class="q-ml-lg","Download",icon="download", @on(:click, :download_event))])
        cell(["Download Info: {{downloadinfo}}"]) ]),@showif("montre_moi_tirer"))
    ]

end

@page("/reactive", ui)
end

# Stipple.range(
        # 0:5:100,
        # :Range_Markers_r,
        # markers = true,
        # var":marker-labels" = "selectionintervalle",
        # color = "secondary"
        # )