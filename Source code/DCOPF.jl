#Project: SOCP Formulation of OPF


#Warning: Paths used here have Linux or Mac format (i.e. /..../), which is different to Windows format ("\\.....\\")
println(string("  "))
println(string("      .................. #  # ...................      "))
println(string("      DCOPF Formulation      "))
println(string("      Author: Paul Serna-Torre, Dinah Shi"))
println(string("      Lab: REAM Lab, University of California San Diego"))
println(string("      Language: Julia 1.6.3"))
println(string("      .................. #  # ...................      "))
println(string("  "))
println(string("  "))


#Packages
using DataFrames, CSV, CPLEX, JuMP, LinearAlgebra, Dates, Plots, VegaLite, Ipopt
println(string("Julia packages were loaded successfully"))

println(string("Reading input files (csv) and creating dataframes and arrays. . ."))

#Configure the fontfamily of the plots
default(; # Plots defaults
fontfamily="Computer modern",
label="" # only explicit legend entries
)

#Formulation
formulation="DCOPF";

#Load directory path
paths=CSV.read("directory_paths.csv", DataFrame);
folder_name=paths[1,2];
input_folder=paths[2,2];
output_folder=paths[3,2];
renewable_folder=paths[4,2];

#Input data folder
path_inputdata=string(folder_name, "/",input_folder,"/");

#Output data folder
path_output=string(folder_name, "/",output_folder, "/");

#Load Configuration.csv
config=CSV.read(string(path_inputdata,"configuration.csv"), DataFrame);
case_name=config[1,2]
start_case=config[2,2];
final_case=config[3,2];
date_format=config[4,2];
resolution=parse(Int,config[5,2]);
price_p_shed=parse(Float64,config[6,2]);
price_q_shed=parse(Float64,config[7,2]);
tolerance_solver=parse(Float64,config[8,2]);
cores_solver=parse(Float64, config[9,2]);
MVAbase=parse(Float64, config[10,2]);
digits_round=parse(Int, config[11,2]);
consider_bess=parse(Int, config[12,2]);
bess_modeling=parse(Int, config[13,2]);
voltage_bound=parse(Float64,config[14,2]);
consider_ctrans=parse(Int, config[15,2]);

#Set main parameters
start_case=Date(start_case,date_format);
final_case=Date(final_case,date_format);
ndays=Dates.value(final_case-start_case)+1;
deltaT=resolution/60;
nper=Int(24/deltaT);

#Load bus.csv
bus=CSV.read(string(path_inputdata,"bus.csv"), DataFrame)[:,1:6];
nbuses=size(bus,1);
insertcols!(bus, 1, :order => 1:nbuses);

#Create a time chart.csv to know the time horizon and duration of each period
time_chart=Array{Any}(nothing,ndays*nper+1,8);
time_chart[1,1:8]=["nday" "nper" "year" "month" "day" "hour" "min" "date"];

    for d in 1:ndays, t in 1:nper    
            time_chart[t+1+(d-1)*nper,1]=d;
            time_chart[t+1+(d-1)*nper,2]=t;
            time_chart[t+1+(d-1)*nper,3]=Dates.year(start_case+Dates.Day(d-1));
            time_chart[t+1+(d-1)*nper,4]=Dates.month(start_case+Dates.Day(d-1));
            time_chart[t+1+(d-1)*nper,5]=Dates.day(start_case+Dates.Day(d-1));
            time_chart[t+1+(d-1)*nper,6]=Dates.hour(DateTime(start_case)+Dates.Minute(resolution*(t-1)));
            time_chart[t+1+(d-1)*nper,7]=Dates.minute(DateTime(start_case)+Dates.Minute(resolution*(t-1)));        
            time_chart[t+1+(d-1)*nper,8]=DateTime(start_case+Dates.Day(d-1))+Dates.Minute(resolution*(t-1));
    end
    time_chart_df = DataFrame(time_chart[2:end,:], :auto);
    rename!(time_chart_df, Symbol.(time_chart[1,:]));
    CSV.write(string(path_output,"time_chart.csv"), time_chart_df);

    time_chart_df[!,:nday] = convert.(Int, time_chart_df[!,:nday]);
    time_chart_df[!,:nper] = convert.(Int, time_chart_df[!,:nper]);
    time_chart_df[!,:year] = convert.(Int, time_chart_df[!,:year]);
    time_chart_df[!,:month] = convert.(Int, time_chart_df[!,:month]);
    time_chart_df[!,:day] = convert.(Int, time_chart_df[!,:day]);
    time_chart_df[!,:hour] = convert.(Int, time_chart_df[!,:hour] );
    time_chart_df[!,:min] = convert.(Int, time_chart_df[!,:min]);

#Upload plant.csv
plant=CSV.read(string(path_inputdata,"plant.csv"), DataFrame)[:,[1,2,5,9,10,11,27,29,30]];
ngen=size(plant,1);

#Non-dispatchable units
dis_units=plant[ ((plant.type .!= "wind") .& (plant.type .!= "solar")) .& (plant.status .==1), :];
ngendis=size(dis_units,1);
insertcols!(dis_units, 1, :order => 1:ngendis);
gencost(g)=dis_units[g, :GenFuelCost]*dis_units[g, :GenIOB]; 

nondis_units=plant[ (plant.type .== "wind") .| (plant.type .== "solar") .& (plant.status .==1), :];
ngennodis=size(nondis_units,1);
insertcols!(nondis_units, 1, :order => 1:ngennodis);

#Set the path to pull out renewable energy data
path_reprofile=string(path_inputdata, renewable_folder);

#Function that finds the maximum renewable power availability of the generator g on day d at period t
function prnw(g,d,t)
    zone= bus[ bus.bus_agg .== nondis_units[g,:bus_agg], :zone_id][1,1];
    type= nondis_units[g,:type];
    profile=CSV.read(string(path_reprofile,"/",type,"_zone",zone,".csv"),DataFrame,header=4)[:,[2,3]];
    profile.time = DateTime.(profile.local_time, "yyyy-mm-dd H:M");
    #Maximum power availability = maximum capacity * capacity factor
    value=nondis_units[g,:Pmax]*profile[ (profile.time .== time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1]),:electricity][1,1];
        if ismissing(value)
                return 0
        else
                return value
        end
end

#Load load.csv
load=CSV.read(string(path_inputdata,"load.csv"), DataFrame);
load.Date = Date.(load.Date, "mm/dd/yyyy")

#Function that finds the active power demand (MW) of the bus n on day d at period t
function pc(n,d,t)  
    global column=0;
    for i in 1:size(names(load),1)
        if names(load)[i]==string(bus[n, :bus_agg])
           global column = i;
           break
        end
   end
   if column !=0 
    value=load[ (load.Date .== Date(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])) .& 
    (load.Time .== Time(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])),column][1,1];    
   else
    value=0
   end

    if ismissing(value)
            return 0
    else
            return value
    end
end 

#Function that calculate the active power demand (MW) of the bus n on day d at period t
function qc(n,d,t)  
   value=pc(n,d,t)*tan(acos(bus[n,:cos]));
    if ismissing(value)
            return 0
    else
            return value
    end
end 

#Some functions to arrange branch data

#Function that find the branch if (m,n) belongs to branches matrix
function map(m,n,branches)
    a=[m n]; 
    found=0;
    rowfound=0;
    Searchi=Int[ a == [branches[i,1] branches[i,2]] for i=1:size(branches,1) ]
    for i in 1:size(Searchi,1)
        if Searchi[i,1] == 1
            found=1
            rowfound=i
        end
    end
    if found==0
        a=[n m];
        Searchi=Int[ a == [branches[i,1] branches[i,2]] for i=1:size(branches,1) ]
        for i in 1:size(Searchi,1)
            if Searchi[i,1] == 1
                found=1
                rowfound=i
            end
        end
    end
    return rowfound
end

#Function that verify if the branch (m,n) exists
function exist(m,n,branches)
    a=[m n]; 
    found=0;
    rowfound=0;
    Searchi=Int[ a == [branches[i,2] branches[i,3]] for i=1:size(branches,1) ]
    for i in 1:size(Searchi,1)
        if Searchi[i,1] == 1
            found=1
            rowfound=i
        end
    end
  
    return found
end

#Load branch.csv
branch=CSV.read(string(path_inputdata,"branch.csv"), DataFrame)[:,1:8];
branch=branch[branch.status.==1,:];
nbranches=size(branch,1);
insertcols!(branch, 1, :order => 1:nbranches);
mapbranch=zeros(nbranches,3);
for br in 1:nbranches
    mapbranch[br,1:3]=[br bus[bus.bus_agg.==branch[br,:from_bus_agg],:order][1] bus[bus.bus_agg.==branch[br,:to_bus_agg],:order][1]];
end
mapbranch=convert.(Int, mapbranch);
insertcols!(branch, 10, :maxcurrent => branch[:, :rateA]/(bus[(bus.bus_agg.==branch[:, :from_bus_agg][1]),:voltage_kV])*bus[(bus.bus_agg.==branch[:, :from_bus_agg][1]),:voltage_kV]/MVAbase);

#Load battery.csv
bess=CSV.read(string(path_inputdata,"battery.csv"), DataFrame);
nbess=size(bess,1);
insertcols!(bess, 1, :order => 1:nbess);

#Main sets of the optimization problem
N = 1:nbuses;
BR = 1:nbranches;
NPER = 1:nper;
TH = 1:ngendis;
RE = 1:ngennodis;
BESS = 1:nbess;

println(string("All input files (csv) were loaded successfully"))
println(string(" "))
println(string("General characteristics of the system under analysis: "))
println(string("Case name: ", case_name))
println(string("Buses: ", nbuses))
println(string("Transmission lines: ", nbranches))
println(string("Resolution (periods per day): ", nper))
println(string("Dispatchable generators: ", ngendis))
println(string("Non-dispatchable generators: ", ngennodis))
println(string("Consider transmission capacity: ", consider_ctrans==1 ? "Yes" : "No"))
println(string("Consider BESS?: ", consider_bess==1 ? "Yes" : "No"))
if consider_bess==1
    println(string("Batteries (BESS): ", nbess))
    println(string("Accurate BESS modeling activated?: ", bess_modeling==1 ? "Yes" : "No"))
end

println(string(" "))

#Admittance matrix Y
modbr=zeros(nbranches);
phibr=zeros(nbranches);
G=zeros(nbuses,nbuses);
B=zeros(nbuses,nbuses);
modY=zeros(nbuses,nbuses);
phiY=zeros(nbuses,nbuses);

for br in BR
    modbr[br]=1/(branch[br,:r]^2+branch[br,:x]^2)^0.5;
    phibr[br]=-atan(branch[br,:x],branch[br,:r]);
end

for i in eachrow(mapbranch)
    l=i[1]; n=i[2]; m=i[3];

    G[n,m]= G[n,m] - branch[l, :status]*modbr[l]*cos(phibr[l]);
    G[m,n]= G[m,n] - branch[l, :status]*modbr[l]*cos(phibr[l]);
    G[n,n]= G[n,n] + branch[l, :status]*modbr[l]*cos(phibr[l]);
    G[m,m]= G[m,m] + branch[l, :status]*modbr[l]*cos(phibr[l]);

    B[n,m]= B[n,m] - branch[l, :status]*modbr[l]*sin(phibr[l]);
    B[m,n]= B[m,n] - branch[l, :status]*modbr[l]*sin(phibr[l]);
    B[n,n]= B[n,n] + branch[l, :status]*modbr[l]*sin(phibr[l]);
    B[m,m]= B[m,m] + branch[l, :status]*modbr[l]*sin(phibr[l]);
    
end

modY[:,:]=[(G[n,m]^2+B[n,m]^2)^0.5 for n in N, m in N];

phiY[:,:]=[ modY[n,m]!=0 ?  atan(B[n,m]/G[n,m]) : 0   for n in N, m in N];

println(string("The admittance matrix was computed successfully "))

d=1;

    println(string("    "))
    println(string("Day ",d))

    #Creat a folder for OPF results of each day
    global daterun=start_case+Dates.Day(d-1);
    global opffolder=string(path_output,"Results ", 
                            string(case_name)," ",
                            string(formulation)," ",                        
                            (consider_bess==0 ? "No BESS" : "BESS")
                            );

    #Create folder

    mkpath(opffolder);

    println(string("Generating the load profile and renewable energy generation..."))

    global pcc=zeros(nbuses,nper);
    global qcc=zeros(nbuses,nper);
    global premax=zeros(ngennodis,nper);

    maxre_results=Array{Any}(nothing,nper+3,ngennodis+1);
    maxre_results[1:3,1]=["plant_id" "bus_agg" "type"]; 

    loadp_results=Array{Any}(nothing,nper+2,nbuses+1);
    loadp_results[1:2,1]=["order" "bus_agg"];

    loadq_results=Array{Any}(nothing,nper+2,nbuses+1);
    loadq_results[1:2,1]=["order" "bus_agg"];

        for t in NPER
            maxre_results[t+3,1]=t;
            for g in RE
                maxre_results[1:3,g+1]=[nondis_units[g, :plant_id][1] nondis_units[g, :bus_agg][1] nondis_units[g, :type][1]] ;
                premax[g,t]=prnw(g,d,t);
                maxre_results[t+3,g+1]=premax[g,t];
            end

            loadp_results[t+2,1]=t;
            loadq_results[t+2,1]=t;
            for n in N
                loadp_results[1:2,n+1]=[bus[n, :order][1] bus[n, :bus_agg][1]];
                loadq_results[1:2,n+1]=[bus[n, :order][1] bus[n, :bus_agg][1]];
                pcc[n,t]=pc(n,d,t);
                qcc[n,t]=qc(n,d,t);
                loadp_results[t+2,n+1]=pcc[n,t];              
                loadq_results[t+2,n+1]=qcc[n,t];     
             end
        end

    maxre_results_df = DataFrame(maxre_results[2:end,:], :auto);
    rename!(maxre_results_df, Symbol.(maxre_results[1,:]));
    CSV.write(string(opffolder,"/","renewable energy availability(MW).csv"), maxre_results_df);
 
    loadp_results_df = DataFrame(loadp_results[2:end,:], :auto);
    rename!(loadp_results_df, Symbol.(loadp_results[1,:]));
    CSV.write(string(opffolder,"/","load(MW).csv"), loadp_results_df);

    loadq_results_df = DataFrame(loadq_results[2:end,:], :auto);
    rename!(loadq_results_df, Symbol.(loadq_results[1,:]));
    CSV.write(string(opffolder,"/","load(MVAR).csv"), loadq_results_df);

    println(string("The files renewable generation availability.csv and load.csv were printed successfully"))

    #Power Flow Model
    global pthh=zeros(ngendis,nper);
    global qthh=zeros(ngendis,nper);
    global pree=zeros(ngennodis,nper);
    global pshedd=zeros(nbuses,nper);
    global qshedd=zeros(nbuses,nper);
    global cc=zeros(nbuses,nbuses,nper);
    global ss=zeros(nbuses,nbuses,nper);
    global thetaa=zeros(nbuses,nper);
    global pdiss=zeros(nbess,nper);
    global pchh=zeros(nbess,nper);
    global ee=zeros(nbess,nper);
    global lmpp=zeros(nbuses,nper);
    global lmpq=zeros(nbuses,nper);
    global dualcurr=zeros(nbranches,nper);
    global ibb=zeros(nbess,nper);
        
    println(string("DC Power Flow"))
    println(string("Day ", d, " Running..."))

    PowerFlow = JuMP.Model(CPLEX.Optimizer);
    @variable(PowerFlow, Cost);
    @variable(PowerFlow, 0 <= pth[TH, NPER]);
    @variable(PowerFlow, 0 <= pre[RE, NPER]);
    @variable(PowerFlow, 0 <= pshed[N, NPER]);
    @variable(PowerFlow, -1.5708<=theta[N, NPER]<=1.5708);


    if consider_bess==1
        @variable(PowerFlow, 0 <= e[BESS, NPER]);
        @variable(PowerFlow, 0 <= pch[BESS, NPER]);
        @variable(PowerFlow, 0 <= pdis[BESS, NPER]);
    end 

    #Objective Function
    @objective(PowerFlow, Min, Cost);

    @constraint(PowerFlow, reference[t in NPER], 
    theta[1,t]==0);

    #Active power balance
    @constraint(PowerFlow, pbalance[n in N, t in NPER], 
            deltaT*(sum( dis_units[g, :bus_agg]==bus[n,:bus_agg] ? pth[g,t] : 0 for g in TH) 
            + sum( nondis_units[g, :bus_agg]==bus[n,:bus_agg] ? pre[g,t] : 0 for g in RE) 
            + (bus[n, :Pd]!=0 ? pshed[n,t] : 0) - pcc[n,t]
            + (consider_bess==1 ? sum( bess[b, :bus_agg]==bus[n,:bus_agg] ? (pdis[b,t] - pch[b,t]) : 0 for b in BESS) : 0))
            == deltaT*(-MVAbase*sum( n!=m ? (B[n,m]*(theta[n,t]-theta[m,t])) : 0 for m in N)));

    #Maximum bound of active power generation of dispatchable power plants
    @constraint(PowerFlow, maxpth[g in TH, t in NPER], 
            pth[g,t] <=  dis_units[g, :Pmax]);
    
    #Minimum bound of active power generation of dispatchable power plants
    @constraint(PowerFlow, minpth[g in TH, t in NPER], 
            pth[g,t] >=  dis_units[g, :Pmin]);

    #Maximum renewable power generation of non-dispatchable power plants
    @constraint(PowerFlow, maxpre[g in RE, t in NPER], 
            pre[g,t] <=   premax[g,t]);

    @constraint(PowerFlow, psheddconst[n in N, t in NPER],
            pshed[n,t] <= pcc[n,t]);


    if  consider_ctrans==1
    @constraint(PowerFlow, current_constraint1[br in BR, t in NPER],  
            MVAbase*modbr[br]*sin(phibr[br])*(theta[mapbranch[br,2],t]-theta[mapbranch[br,3],t]) <= branch[br, :rateA])

    @constraint(PowerFlow, current_constraint2[br in BR, t in NPER],
            MVAbase*modbr[br]*sin(phibr[br])*(theta[mapbranch[br,2],t]-theta[mapbranch[br,3],t]) >= -branch[br, :rateA])                        
    end


    if consider_bess==1

        @constraint(PowerFlow, BessState[b in BESS, t in NPER],
            e[b,t]== (t==1 ? (bess[b, :Initial_state]*bess[b, :Emax]) : e[b,t-1])
                     +(pch[b,t]*bess[b, :ChargEff] - pdis[b,t]/bess[b, :DischEff])*deltaT);
    
        @constraint(PowerFlow, BessFin[b in BESS], 
                     e[b,nper]==bess[b, :Final_state]*bess[b, :Emax]);
            
        @constraint(PowerFlow, BessDmax[b in BESS, t in NPER],
                     pdis[b,t] <= bess[b, :Pmax]); 
    
        @constraint(PowerFlow, BessCmax[b in BESS, t in NPER],
                     pch[b,t] <= bess[b, :Pmax]);   
    
                if bess_modeling==0
                    @constraint(PowerFlow, BessEmax[b in BESS, t in NPER],
                    e[b,t] <= bess[b, :Emax]);

                else
                    @constraint(PowerFlow, BessEmax[b in BESS, t in NPER],
                    e[b,t] <= bess[b, :Emax]);
                end
    end

    #Objective function
    @constraint(PowerFlow, z, 
                        Cost==
                        price_p_shed*deltaT*sum(bus[n, :Pd]!=0 ? pshed[n,t] : 0 for n in N, t in NPER) +
                        deltaT*sum(gencost(g)*pth[g,t] for g in TH,t in NPER) 
   #                    ((consider_bess==1) ? (10^-2*sum(pdis[b,t]/MVAbase for b in BESS, t in NPER)) : 0) 
                        );

    JuMP.optimize!(PowerFlow);
    global status1=string(primal_status(PowerFlow));
    global status2=string(termination_status(PowerFlow));
    
    println(string("Day ", d," - ", primal_status(PowerFlow)) ,"-", termination_status(PowerFlow));
    
    if status1=="FEASIBLE_POINT"  
        global  pthh=JuMP.value.(pth);
        global  pree=JuMP.value.(pre);
        global  pshedd=JuMP.value.(pshed);
        global  thetaa=JuMP.value.(theta);
        global  lmpp=JuMP.dual.(pbalance);

        if consider_ctrans==1
            global  dualcurr=JuMP.dual.(current_constraint1).+JuMP.dual.(current_constraint2);
        end
        
        if consider_bess==1
         global  pdiss=round.(JuMP.value.(pdis),digits=digits_round);
         global  pchh=round.(JuMP.value.(pch), digits=digits_round);
         global  ee=round.(JuMP.value.(e), digits=digits_round);    
        end

        cc[:,:,:].=1;
        for t in NPER
            for n in N
                for m in N
                    if B[n,m]!=0
                    ss[n,m,t]=thetaa[n,t]-thetaa[m,t];
                    # it should be thetaa[m,t]-thetaa[n,t] (I don't understand);
                    end
                end
            end     
        end
    end

     #Call the file Print_CSVfiles.jl and Print_graphs.jl
     println(string(" "))
     println(string("Creating output files . . ."))

     include(string(folder_name,"/Source code/","Print_CSVfiles.jl"))     

     include(string(folder_name,"/Source code/","Print_graphs.jl"))     




