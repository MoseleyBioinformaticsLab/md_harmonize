Search.setIndex({docnames:["api","guide","index","license","tutorial"],envversion:{"sphinx.domains.c":1,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":1,"sphinx.domains.index":1,"sphinx.domains.javascript":1,"sphinx.domains.math":2,"sphinx.domains.python":1,"sphinx.domains.rst":1,"sphinx.domains.std":1,"sphinx.ext.intersphinx":1,"sphinx.ext.todo":2,"sphinx.ext.viewcode":1,sphinx:56},filenames:["api.rst","guide.rst","index.rst","license.rst","tutorial.rst"],objects:{"":{md_harmonize:[0,0,0,"-"]},"md_harmonize.KEGG_database_scraper":{curate_molfile:[0,1,1,""],download:[0,1,1,""],entry_list:[0,1,1,""],update_entity:[0,1,1,""]},"md_harmonize.KEGG_parser":{RpairParser:[0,2,1,""],compound_pair_mappings:[0,1,1,""],create_atom_mappings:[0,1,1,""],create_compound_kcf:[0,1,1,""],create_reactions:[0,1,1,""],kegg_data_parser:[0,1,1,""],kegg_kcf_parser:[0,1,1,""],parse_equation:[0,1,1,""],reaction_center:[0,2,1,""]},"md_harmonize.KEGG_parser.RpairParser":{__init__:[0,3,1,""],combine_atom_mappings:[0,3,1,""],construct_component:[0,3,1,""],count_changed_atom_identifiers:[0,3,1,""],create_reaction_centers:[0,3,1,""],detect_components:[0,3,1,""],find_center_atoms:[0,3,1,""],find_target_atom:[0,3,1,""],generate_atom_mappings:[0,3,1,""],generate_kat_neighbors:[0,3,1,""],get_center_list:[0,3,1,""],map_atom_by_colors:[0,3,1,""],map_components:[0,3,1,""],map_whole_compound:[0,3,1,""],pair_components:[0,3,1,""],preliminary_atom_mappings_check:[0,3,1,""],remove_different_bonds:[0,3,1,""],validate_component_atom_mappings:[0,3,1,""]},"md_harmonize.KEGG_parser.reaction_center":{__getnewargs__:[0,3,1,""],__new__:[0,3,1,""],__repr__:[0,3,1,""],difference:[0,4,1,""],i:[0,4,1,""],kat:[0,4,1,""],label:[0,4,1,""],match:[0,4,1,""]},"md_harmonize.MetaCyc_parser":{atom_mappings_parser:[0,1,1,""],create_reactions:[0,1,1,""],generate_one_to_one_mappings:[0,1,1,""],reaction_parser:[0,1,1,""],reaction_side_parser:[0,1,1,""]},"md_harmonize.aromatics":{AromaticManager:[0,2,1,""]},"md_harmonize.aromatics.AromaticManager":{__init__:[0,3,1,""],add_aromatic_substructures:[0,3,1,""],construct_aromatic_entity:[0,3,1,""],decode:[0,3,1,""],detect_aromatic_substructures:[0,3,1,""],detect_aromatic_substructures_timeout:[0,3,1,""],encode:[0,3,1,""],extract_aromatic_substructures:[0,3,1,""],fuse_cycles:[0,3,1,""],indigo_aromatic_bonds:[0,3,1,""],indigo_aromatize:[0,3,1,""],kegg_aromatize:[0,3,1,""]},"md_harmonize.compound":{Atom:[0,2,1,""],Bond:[0,2,1,""],Compound:[0,2,1,""]},"md_harmonize.compound.Atom":{__init__:[0,3,1,""],add_neighbors:[0,3,1,""],clone:[0,3,1,""],color_atom:[0,3,1,""],remove_neighbors:[0,3,1,""],reset_color:[0,3,1,""],update_atom_number:[0,3,1,""],update_cycle:[0,3,1,""],update_kat:[0,3,1,""],update_stereochemistry:[0,3,1,""],update_symbol:[0,3,1,""]},"md_harmonize.compound.Bond":{__init__:[0,3,1,""],clone:[0,3,1,""],update_bond_type:[0,3,1,""],update_first_atom:[0,3,1,""],update_second_atom:[0,3,1,""],update_stereochemistry:[0,3,1,""]},"md_harmonize.compound.Compound":{__init__:[0,3,1,""],backbone_color_identifier:[0,3,1,""],break_cycle:[0,3,1,""],calculate_bond_stereochemistry:[0,3,1,""],calculate_distance_to_r_groups:[0,3,1,""],calculate_y_coordinate:[0,3,1,""],circular_pair_relationship:[0,3,1,""],circular_pair_relationship_helper:[0,3,1,""],collect_atomic_weights_of_neighbors:[0,3,1,""],color_compound:[0,3,1,""],color_groups:[0,3,1,""],color_h:[0,3,1,""],color_metal:[0,3,1,""],compare_branch_weights:[0,3,1,""],compare_chemical_details:[0,3,1,""],compare_chemical_details_with_mapping:[0,3,1,""],composition:[0,3,1,""],connected_components:[0,3,1,""],contains_r_groups:[0,3,1,""],curate_invalid_n:[0,3,1,""],curate_invalid_symmetric_atoms:[0,3,1,""],define_bond_stereochemistry:[0,3,1,""],detect_abnormal_atom:[0,3,1,""],determine_relationship:[0,3,1,""],distance_matrix:[0,3,1,""],encode:[0,3,1,""],extract_aromatic_bonds:[0,3,1,""],extract_double_bond_connecting_cycle:[0,3,1,""],find_critical_atom_in_cycle:[0,3,1,""],find_cycles:[0,3,1,""],find_cycles_helper:[0,3,1,""],find_double_bond_linked_atom:[0,3,1,""],find_mappings:[0,3,1,""],find_mappings_reversed:[0,3,1,""],first_round_color:[0,3,1,""],formula:[0,3,1,""],generate_atom_color_with_neighbors:[0,3,1,""],generate_atom_mapping_by_atom_color:[0,3,1,""],generate_atom_zero_layer_color:[0,3,1,""],get_chemical_details:[0,3,1,""],get_next_layer_neighbors:[0,3,1,""],h_color_identifier:[0,3,1,""],h_index:[0,3,1,""],has_isolated_atoms:[0,3,1,""],heavy_atoms:[0,3,1,""],index_of_heavy_atoms:[0,3,1,""],invalid_symmetric_atoms:[0,3,1,""],map_r_correspondents:[0,3,1,""],map_resonance:[0,3,1,""],map_resonance_helper:[0,3,1,""],metal_color_identifier:[0,3,1,""],metal_index:[0,3,1,""],molfile_name:[0,3,1,""],name:[0,3,1,""],optimal_mapping_with_r:[0,3,1,""],optimal_resonant_mapping:[0,3,1,""],r_groups:[0,3,1,""],reset_color:[0,3,1,""],restore_cycle:[0,3,1,""],same_structure_relationship:[0,3,1,""],separate_connected_components:[0,3,1,""],structure_matrix:[0,3,1,""],update_aromatic_bond_type:[0,3,1,""],update_atom_symbol:[0,3,1,""],update_color_tuple:[0,3,1,""],validate_mapping_with_r:[0,3,1,""],with_r_pair_relationship:[0,3,1,""],with_r_pair_relationship_helper:[0,3,1,""]},"md_harmonize.harmonization":{CompoundHarmonizationManager:[0,2,1,""],HarmonizationManager:[0,2,1,""],HarmonizedCompoundEdge:[0,2,1,""],HarmonizedEdge:[0,2,1,""],HarmonizedReactionEdge:[0,2,1,""],ReactionHarmonizationManager:[0,2,1,""],harmonize_compound_list:[0,1,1,""],harmonize_reaction_list:[0,1,1,""]},"md_harmonize.harmonization.CompoundHarmonizationManager":{__init__:[0,3,1,""],add_edge:[0,3,1,""],add_invalid:[0,3,1,""],create_manager:[0,3,1,""],find_compound:[0,3,1,""],get_edge_list:[0,3,1,""],has_visited:[0,3,1,""],remove_edge:[0,3,1,""]},"md_harmonize.harmonization.HarmonizationManager":{__init__:[0,3,1,""],add_edge:[0,3,1,""],create_key:[0,3,1,""],remove_edge:[0,3,1,""],save_manager:[0,3,1,""],search:[0,3,1,""]},"md_harmonize.harmonization.HarmonizedCompoundEdge":{__init__:[0,3,1,""],pair_atom_mappings:[0,3,1,""],reversed_mappings:[0,3,1,""]},"md_harmonize.harmonization.HarmonizedEdge":{__init__:[0,3,1,""],pair_relationship:[0,3,1,""],reversed_relationship:[0,3,1,""]},"md_harmonize.harmonization.HarmonizedReactionEdge":{__init__:[0,3,1,""]},"md_harmonize.harmonization.ReactionHarmonizationManager":{__init__:[0,3,1,""],compare_ecs:[0,3,1,""],compound_mappings:[0,3,1,""],determine_relationship:[0,3,1,""],harmonize_reaction:[0,3,1,""],jaccard:[0,3,1,""],match_unmapped_compounds:[0,3,1,""],one_to_one_compound_mappings:[0,3,1,""],unmapped_compounds:[0,3,1,""]},"md_harmonize.reaction":{Reaction:[0,2,1,""]},"md_harmonize.reaction.Reaction":{__init__:[0,3,1,""],name:[0,3,1,""]},md_harmonize:{KEGG_database_scraper:[0,0,0,"-"],KEGG_parser:[0,0,0,"-"],MetaCyc_parser:[0,0,0,"-"],aromatics:[0,0,0,"-"],compound:[0,0,0,"-"],harmonization:[0,0,0,"-"],reaction:[0,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","function","Python function"],"2":["py","class","Python class"],"3":["py","method","Python method"],"4":["py","attribute","Python attribute"]},objtypes:{"0":"py:module","1":"py:function","2":"py:class","3":"py:method","4":"py:attribute"},terms:{"break":0,"case":0,"class":0,"final":0,"float":0,"function":[0,1,4],"int":0,"new":0,"public":[1,4],"return":0,"short":0,"static":0,"true":0,"try":0,AND:3,ARE:3,And:0,BUT:3,FOR:3,For:0,NOT:3,One:0,SUCH:3,THE:3,The:[1,2],Then:0,There:0,USE:3,Use:1,Used:0,__getnewargs__:0,__init__:0,__new__:0,__repr__:0,_cl:0,about:[0,1],abov:[0,3],access:[0,1],account:0,acetyl:0,acetylcitrullin:0,acetylglutam:0,acetyltransferas:0,achiev:[0,1],acid:0,across:[0,1],add:[0,4],add_aromatic_substructur:0,add_edg:0,add_invalid:0,add_neighbor:0,added:[0,4],adding:0,addit:0,advis:3,after:0,aldehyd:0,aldol:0,algorithm:0,alia:0,all:[0,3],alreadi:0,also:0,amino:0,ani:[0,3],api:[1,2],appli:0,arginin:0,argininosuccin:0,aris:3,aromat:[1,2],aromatic_cycl:0,aromatic_manag:[1,4],aromatic_manager_fil:[1,4],aromatic_structur:0,aromatic_substructur:0,aromaticmanag:0,associ:0,assum:0,atom:[0,1,4],atom_atom_mapping_numb:0,atom_c:0,atom_forming_double_bond:0,atom_in_cycl:0,atom_index:0,atom_map:0,atom_mapping_fil:0,atom_mapping_text:0,atom_mappings_pars:0,atom_numb:0,atom_o:0,atom_oo:0,atom_stereo:0,atom_stereo_par:0,atom_symbol:0,atoms_to_color:0,attent:0,avail:[1,4],avoid:0,backbon:0,backbone_color_identifi:0,balanc:0,base:[0,4],basic:[0,2],bass:[0,4],been:0,belong:0,below:3,between:[0,1,4],bifunct:0,binari:3,biosynthesi:0,bond:0,bond_stereo:0,bond_topolog:0,bond_typ:0,bool:0,both:0,branch:0,breadth:0,break_cycl:0,broken:0,bsd:3,bunch:0,busi:3,c00001:0,c00003:0,c00004:0,c00010:0,c00010_c00024:0,c00013:0,c00024:0,c00025:0,c00025_c00624:0,c00029:0,c00080:0,c00167:0,c00624:0,c1c:0,c8x:0,c8y:0,calcul:[0,1],calculate_bond_stereochemistri:0,calculate_distance_to_r_group:0,calculate_y_coordin:0,call:0,can:[0,1,4],candid:0,cannot:0,captur:0,care:0,caspi:0,caus:[0,3],cco:0,center:0,center_atom_index:0,center_atom_numb:0,chain:0,chang:0,character:0,charg:0,check:0,chemic:0,chemistri:0,choos:0,circular:0,circular_pair_relationship:0,circular_pair_relationship_help:0,citat:3,clean:0,clear:3,clone:[0,1],coa:0,code:[2,3],coeffici:0,collect:0,collect_atomic_weights_of_neighbor:0,color:0,color_0:0,color_atom:0,color_compound:0,color_group:0,color_h:0,color_met:0,com:1,combin:0,combine_atom_map:0,comma:4,command:1,commiss:0,compar:0,compare_branch_weight:0,compare_chemical_detail:0,compare_chemical_details_with_map:0,compare_ec:0,comparison:0,compart:0,compon:0,component_atom_map:0,composit:0,compound:[1,2],compound_dict:0,compound_dict_list:0,compound_harmonization_manag:0,compound_map:0,compound_nam:0,compound_pair:0,compound_pair_map:0,compoundharmonizationmanag:0,concaten:0,condit:3,confus:0,connect:0,connected_compon:0,consequenti:3,consid:0,construct:[0,1],construct_aromatic_ent:0,construct_compon:0,contain:0,contains_r_group:0,content:2,contract:3,contributor:3,coordin:0,copi:[0,4],copyright:3,core:0,correspond:0,count:0,count_changed_atom_identifi:0,counterpart:0,cpd:0,creat:[0,1],create_atom_map:0,create_compound_kcf:0,create_kei:0,create_manag:0,create_react:0,create_reaction_cent:0,credit:0,criteria:0,critic:0,critical_atom:0,ctfile:[0,1],cur_layer_neighbor:0,curat:0,curate_invalid_n:0,curate_invalid_symmetric_atom:0,curate_molfil:0,current:0,cutoff:0,cycl:0,cycle_statu:0,cython:[0,1],damag:3,data:[0,2,3],databas:[0,1],database_nam:[1,4],dblink:0,decod:0,defin:0,define_bond_stereochemistri:0,definit:[0,4],depend:2,depth:0,deriv:[0,3],describ:0,descript:[0,2],design:0,detail:[0,1,4],detect:[0,1,4],detect_abnormal_atom:0,detect_aromatic_substructur:0,detect_aromatic_substructures_timeout:0,detect_compon:0,determin:0,determine_relationship:0,dict:0,dictionari:0,differ:0,dijkstra:0,direct:[0,3],directli:[0,1],directori:[0,4],disassembl:0,disclaim:3,discov:4,disjoint:0,distanc:0,distance_matrix:0,distinguish:0,distribut:3,divid:0,docopt:1,document:[3,4],doe:0,doi:[1,4],don:0,doubl:0,doubli:0,down:0,download:[0,1],due:0,duplic:0,each:[0,4],ecs:0,edg:0,edge_typ:0,either:[0,4],element:0,els:0,empti:0,encod:0,end:0,endors:3,ensur:0,entiti:0,entri:0,entry_list:0,enumer:0,environ:0,enzym:0,enzymat:0,equal:0,equat:0,equival:0,error:0,evalu:0,even:3,event:3,everi:0,exact_change_flag:0,exampl:0,exclud:0,excluded_index:0,execut:0,exemplari:3,exist:0,expand:1,explain:4,express:3,extra:0,extract:0,extract_aromatic_bond:0,extract_aromatic_substructur:0,extract_double_bond_connecting_cycl:0,extrem:0,facilit:0,fals:0,fetch:0,fewer:0,field:0,file:[0,1],file_path:0,filenam:0,find:0,find_center_atom:0,find_compound:0,find_critical_atom_in_cycl:0,find_cycl:0,find_cycles_help:0,find_double_bond_linked_atom:0,find_map:0,find_mappings_revers:0,find_target_atom:0,first:[0,4],first_atom_numb:0,first_round_color:0,fit:3,floyd:0,follow:[0,1,3,4],form:[0,1,3],format:0,formula:0,from:[0,1,3],from_index:0,from_sid:0,full:0,fuse:0,fuse_cycl:0,gener:[0,4],generate_atom_color_with_neighbor:0,generate_atom_map:0,generate_atom_mapping_by_atom_color:0,generate_atom_zero_layer_color:0,generate_kat_neighbor:0,generate_one_to_one_map:0,geometr:0,get:[0,2],get_center_list:0,get_chemical_detail:0,get_edge_list:0,get_next_layer_neighbor:0,git:1,github:1,given:0,glutam:0,good:3,grant:3,graph:0,group:0,guid:2,h0design:0,h_color_identifi:0,h_index:0,half:0,haonitro:0,harmon:[1,2],harmonizationmanag:0,harmonizationmang:0,harmonize_compound:[1,4],harmonize_compound_list:0,harmonize_react:[0,1,4],harmonize_reaction_list:0,harmonized_edg:0,harmonizedcompoundedg:0,harmonizededg:0,harmonizedreactionedg:0,has:[0,1],has_isolated_atom:0,has_visit:0,have:0,heavi:0,heavier:0,heavy_atom:0,heavy_sid:0,help:1,here:0,hierarchi:4,hmdb:[1,4],holder:3,howev:[0,3],http:[1,4],huan:3,hunter:3,hydrogen:0,hydrogen_count:0,hydroxi:0,hydroxylamin:0,idea:0,identifi:0,idx:0,ignor:0,implement:0,impli:3,in_cycl:0,incident:3,includ:[0,3,4],incorpor:0,index:[0,2],index_of_heavy_atom:0,indic:[0,4],indigo:[0,1,4],indigo_aromat:0,indigo_aromatic_bond:0,indirect:3,inform:0,initi:0,initialize_compound:[1,4],initialize_react:[1,4],instal:2,instanc:0,integ:0,intend:4,intercept:0,interchang:0,interfac:1,interrupt:3,invalid:0,invalid_symmetric_atom:0,inversion_retention_flag:0,involv:0,isol:0,isotop:0,isotope_resolv:0,issu:0,ith:0,its:[0,3],jaccard:0,jin:3,json:1,jsonpickl:[0,1],just:0,k00618:0,k00619:0,k00620:0,k11067:0,k14681:0,k14682:0,k22476:0,k22477:0,k22478:0,kat:0,kcf:[0,4],kcf_cpd:0,kcf_file:0,keep:4,kegg:[0,1],kegg_aromat:0,kegg_data_pars:0,kegg_database_scrap:2,kegg_kcf_pars:0,kegg_pars:2,kei:0,kinas:0,label:0,layer:0,lead:0,least:0,left:0,left_cent:0,left_compon:0,left_removed_bond:0,length:0,level:0,liabil:3,liabl:3,librari:1,licens:2,lie:0,light:0,light_sid:0,like:0,limit:[0,3],line:[0,1],linear:0,link:0,linkag:0,list:[0,3],local:[0,4],loos:0,loss:3,lyas:0,m00028:0,m00845:0,mai:3,main:0,mainli:0,major:0,make:0,manag:[0,1,4],mani:0,map:[0,1,4],map_atom_by_color:0,map_compon:0,map_r_correspond:0,map_reson:0,map_resonance_help:0,map_whole_compound:0,mass:0,mass_differ:0,match:0,match_unmapped_compound:0,materi:3,matrix:0,max:0,md_harmon:[1,2],mean:0,meet:0,mention:0,merchant:3,met:3,metabol:[0,1,4],metabolit:0,metacyc:[0,1,4],metacyc_pars:2,metal:0,metal_color_identifi:0,metal_index:0,method:0,minim:0,modif:3,modifi:3,modul:[0,2],mol:0,molecul:0,molfil:[0,1],molfile_nam:0,more:0,moselei:3,moseleybioinformaticslab:1,most:0,multipl:[0,4],multiprocess:1,must:3,n5y:0,name:[0,3,4],name_1:0,name_2:0,namedtupl:0,ndarrai:0,nearest:0,need:0,neglig:3,neighbor:0,neither:3,newli:[0,4],next:0,nice:0,nitrit:0,none:0,nor:3,normal:0,note:[0,4],notic:3,nth:0,number:0,numpi:[0,1],o1c:0,o2c:0,object:1,obsolet:0,occur:0,occurr:0,often:0,onc:0,one:0,one_chemical_detail:0,one_compound:0,one_ec:0,one_r:0,one_react:0,one_sid:0,one_side_left:0,one_to_one_compound_map:0,one_to_one_map:0,onli:0,optim:[0,1],optimal_mapping_with_r:0,optimal_resonant_map:0,option:[0,1],order:0,org:[1,4],ornithin:0,orphan:0,ortholog:0,other:[0,1,3],other_compound:0,other_ec:0,other_react:0,other_sid:0,other_side_left:0,otherwis:[1,3],out:3,output:1,outsid:0,oxocarboxyl:0,oxygen:0,p1b:0,packag:[0,1],page:2,pai:0,pair:0,pair_atom_map:0,pair_compon:0,pair_relationship:0,paramet:0,pari:0,pariti:0,pars:[0,1,4],parse_equ:0,parse_kegg_atom:[1,4],part:0,parti:3,particular:3,patent:3,path:0,pathwai:0,pattern:0,pebbl:1,period:0,permiss:3,permit:3,physiolog:0,pickl:[0,1],piec:0,pip:1,plain:0,plane:0,pleas:[0,4],posit:0,possibl:[0,3],potenti:0,pre:[1,4],preliminary_atom_mappings_check:0,prior:3,prioriti:0,process:0,procur:3,product:3,profit:3,promot:3,proper:3,properti:0,proton:0,provid:[0,1,3,4],publish:3,purpos:3,python3:1,python:1,r00259:0,r_distanc:0,r_group:0,raw:[1,4],rc00004:0,rc00064:0,rclass:[0,4],rclass_definit:0,rclass_directori:0,rclass_nam:0,rdm:0,react:0,reactant:0,reacting_center_statu:0,reaction:[1,2],reaction_cent:0,reaction_directori:0,reaction_fil:0,reaction_list:0,reaction_nam:0,reaction_pars:0,reaction_sid:0,reaction_side_pars:0,reaction_text:0,reactionharmonizationmanag:0,recolor:0,redistribut:3,redox:0,reduct:0,redund:0,refer:2,region:0,relationship:0,relev:0,remov:0,remove_different_bond:0,remove_edg:0,remove_neighbor:0,removed_bond:0,repositori:1,repres:0,represent:[0,1],reproduc:3,requir:[0,1],reserv:3,reset:0,reset_color:0,reson:0,respons:0,rest:0,restor:0,restore_cycl:0,result:1,retain:3,revers:0,reversed_map:0,reversed_relationship:0,rhea:0,right:[0,3],right_cent:0,right_compon:0,right_removed_bond:0,ring:0,rn00220:0,rn01100:0,rn01110:0,rn01210:0,rn01230:0,roughli:0,round:0,rpairpars:0,run:1,rxn:0,s2x:0,same:0,same_structure_relationship:0,save:[0,1],save_fil:[1,4],save_manag:0,screen:1,search:[0,2],second:0,second_atom_numb:0,secondari:0,self:0,separ:[0,4],separate_connected_compon:0,serializ:1,servic:3,set:0,sever:0,shall:3,share:0,short_circuit:0,should:0,show:1,side:0,simpl:0,simpli:0,sinc:0,singl:0,slope:0,softwar:3,solv:0,some:0,sourc:[0,2,3],special:3,specif:[0,3],speed:1,split:1,sri:0,stai:0,standard:[1,4],start:0,statu:0,std:0,step:0,stereo:0,stereo_care_box:0,stereochemistri:0,stop:0,store:0,str:0,strict:3,string:0,structur:0,structure_matrix:0,sub_directori:0,subdirectori:0,subject:3,subset:0,substitut:3,substructur:0,successfulli:0,suffix:0,sugar:0,suggest:0,support:1,suppos:0,sure:0,symbol:0,symmetr:0,synthas:0,tabl:0,take:0,taken:0,target:0,target_url:0,test:0,text:0,than:0,the_oth:0,the_other_chemical_detail:0,the_other_compound:0,the_other_differ:0,the_other_match:0,thei:0,them:0,theori:3,therefor:0,thi:[0,1,3,4],thing:0,this_compound:0,those:0,three:0,through:0,time:0,timeout:0,to_index:0,to_sid:0,togeth:0,tool:1,topolog:0,tort:3,total:0,transform:0,tripl:0,tupl:0,tutori:[1,2],two:0,type:0,under:1,union:0,uniqu:0,unmap:0,unmapped_compound:0,unmapped_count:0,unmatch:0,until:0,updat:[0,4],update_aromatic_bond_typ:0,update_atom_numb:0,update_atom_symbol:0,update_bond_typ:0,update_color_tupl:0,update_cycl:0,update_ent:0,update_first_atom:0,update_kat:0,update_second_atom:0,update_stereochemistri:0,update_symbol:0,updated_symbol:0,url:[0,1],usag:[2,4],use:[0,1,3],used:[0,3,4],user:2,using:[0,4],valenc:0,valid:0,validate_component_atom_map:0,validate_mapping_with_r:0,variabl:0,veri:0,version:1,via:0,visit:0,wai:3,warranti:3,warshal:0,water:0,weight:0,well:0,when:0,where:0,whether:[0,3],which:0,with_r_pair_relationship:0,with_r_pair_relationship_help:0,without:[0,3],word:0,work:3,working_directori:[1,4],written:3,you:1,zenodo:[1,4],zero:0,zero_core_color:0,zero_neighbor_color:0},titles:["The md_harmonize API Reference","User Guide","Welcome to MDH\u2019s documentation!","License","The md_harmonize Tutorial"],titleterms:{The:[0,4],Using:4,across:4,api:[0,4],aromat:[0,4],basic:1,code:1,compound:[0,4],construct:4,curat:4,data:[1,4],databas:4,depend:1,descript:1,document:2,download:4,entiti:4,extract:4,from:4,get:1,guid:1,harmon:[0,4],indic:2,instal:1,kegg:4,kegg_database_scrap:0,kegg_pars:0,licens:3,linux:1,mac:1,md_harmon:[0,4],mdh:2,metacyc_pars:0,molfil:4,option:4,prepar:4,reaction:[0,4],refer:0,represent:4,sourc:1,substructur:4,tabl:2,tutori:4,usag:1,user:1,welcom:2}})