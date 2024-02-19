        #File for printing csv files

        #Results of thermal power
        pth_results=Array{Any}(nothing,nper+1,ngendis+1);
        pth_results[1,1]=string("periods");
        for g in TH, t in NPER
                pth_results[1,g+1]=string("gen ", dis_units[g, :plant_id]);
                pth_results[t+1,[1 g+1]]=[t pthh[g,t]];
        end
        pth_results_df = DataFrame(round.(pth_results[2:end,:], digits=digits_round), :auto);
        rename!( pth_results_df, Symbol.(pth_results[1,:]));
        CSV.write(string(opffolder,"/","thermal power generation(MW).csv"), pth_results_df);

        qth_results=Array{Any}(nothing,nper+1,ngendis+1);
        qth_results[1,1]=string("periods");
        for g in TH, t in NPER
                qth_results[1,g+1]=string("gen ", dis_units[g, :plant_id]);
                qth_results[t+1,[1 g+1]]=[t qthh[g,t]];
        end
        qth_results_df = DataFrame(round.(qth_results[2:end,:], digits=digits_round), :auto);
        rename!( qth_results_df, Symbol.(qth_results[1,:]));
        CSV.write(string(opffolder,"/","thermal power generation(MVAR).csv"), qth_results_df);

        #Results of bus voltage
        v_results=Array{Any}(nothing,nper+1,nbuses+1);
        v_results[1,1]=string("periods");
        for n in N, t in NPER
                v_results[1,n+1]=string("bus ", n);
                v_results[t+1,[1 n+1]]=[t cc[n,n,t]^0.5];
        end
        v_results_df = DataFrame(round.(v_results[2:end,:], digits=digits_round), :auto);
        rename!(v_results_df, Symbol.(v_results[1,:]));
        CSV.write(string(opffolder,"/","voltage(pu).csv"), v_results_df);

        theta_results=Array{Any}(nothing,nper+1,nbuses+1);
        theta_results[1,1]=string("periods");
        for n in N, t in NPER
                theta_results[1,n+1]=string("bus ", n);
                theta_results[t+1,[1 n+1]]=[t 180/3.1416*thetaa[n,t]];
        end
        theta_results_df = DataFrame(round.(theta_results[2:end,:], digits=digits_round), :auto);
        rename!(theta_results_df, Symbol.(theta_results[1,:]));
        CSV.write(string(opffolder,"/","theta(degree).csv"), theta_results_df);

        #Results of p shedding 
        pshed_results=Array{Any}(nothing,nper+1,nbuses+1);
        pshed_results[1,1]=string("periods");
        for n in N, t in NPER
            pshed_results[1,n+1]=string("bus ",n);
            pshed_results[t+1,[1 n+1]]=[t pshedd[n,t]];
        end
        pshed_results_df = DataFrame(round.(pshed_results[2:end,:], digits=digits_round), :auto);
        rename!(pshed_results_df, Symbol.(pshed_results[1,:]));
        CSV.write(string(opffolder,"/","demand response(MW).csv"), pshed_results_df);

        #Results of q shedding 
        qshed_results=Array{Any}(nothing,nper+1,nbuses+1);
        qshed_results[1,1]=string("periods");
        for n in N, t in NPER
            qshed_results[1,n+1]=string("bus ",n);
            qshed_results[t+1,[1 n+1]]=[t qshedd[n,t]];
        end
        qshed_results_df = DataFrame(round.(qshed_results[2:end,:], digits=digits_round), :auto);
        rename!(qshed_results_df, Symbol.(qshed_results[1,:]));
        CSV.write(string(opffolder,"/","demand response(MVAR).csv"), qshed_results_df);

        #Results of renewable power
        if  ngennodis!=0
        pre_results=Array{Any}(nothing,nper+1,ngennodis+1);
        pre_results[1,1]=string("periods");
        for g in RE, t in NPER
                pre_results[1,g+1]=string("gen ", nondis_units[g, :plant_id]);
                pre_results[t+1,[1 g+1]]=[t pree[g,t]];
        end
        pre_results_df = DataFrame(round.(pre_results[2:end,:], digits=digits_round), :auto);
        rename!( pre_results_df, Symbol.(pre_results[1,:]));
        CSV.write(string(opffolder,"/","renewable power generation(MW).csv"), pre_results_df);
       
        
        #Results of curtailment
        curt_results=Array{Any}(nothing,nper+1,ngennodis+1);
        curt_results[1,1]=string("periods");
        for g in RE, t in NPER
                 curt_results[1,g+1]=string("gen ", nondis_units[g, :plant_id]);
                 curt_results[t+1,[1 g+1]]=[t premax[g,t]-pree[g,t]];
        end
        curt_results_df = DataFrame(round.(curt_results[2:end,:], digits=digits_round), :auto);
        rename!( curt_results_df, Symbol.(curt_results[1,:]));
        CSV.write(string(opffolder,"/","curtailment(MW).csv"), curt_results_df );
        end

        #Active power flow across a branch
        pflow_results=Array{Any}(nothing,nper+1,2*nbranches+1);
        pflow_results[1,1]=string("periods");
        plosses=zeros(nbranches,nper)
        for br in BR
            n=mapbranch[br,2];
            m=mapbranch[br,3];    
            pflow_results[1,2*br]=string("id",branch[br,:branch_id]," ",n,"->",m);
            pflow_results[1,2*br+1]=string("id",branch[br,:branch_id]," ",m,"->",n);
            for t in NPER
                pflow_results[t+1,1]=t;
                pflow_results[t+1,2*br]=MVAbase*(modbr[br]*cos(phibr[br])*(cc[n,n,t]-cc[n,m,t])+modbr[br]*sin(phibr[br])*ss[n,m,t]);
                pflow_results[t+1,2*br+1]=MVAbase*(modbr[br]*cos(phibr[br])*(cc[m,m,t]-cc[m,n,t])+modbr[br]*sin(phibr[br])*ss[m,n,t]);
                plosses[br,t]=pflow_results[t+1,2*br]+pflow_results[t+1,2*br+1];
            end
        end
        pflow_results_df = DataFrame(round.(pflow_results[2:end,:], digits=digits_round), :auto);
        rename!(pflow_results_df, Symbol.(pflow_results[1,:]));
        CSV.write(string(opffolder,"/","power flow branches(MW).csv"), pflow_results_df);

        #Reactive power flow across a branch
        qflow_results=Array{Any}(nothing,nper+1,2*nbranches+1);
        qflow_results[1,1]=string("periods");
        qlosses=zeros(nbranches,nper)
        for br in BR
            n=mapbranch[br,2];
            m=mapbranch[br,3];    
            qflow_results[1,2*br]=string("id",branch[br,:branch_id]," ",n,"->",m);
            qflow_results[1,2*br+1]=string("id",branch[br,:branch_id]," ",m,"->",n);
            for t in NPER
                qflow_results[t+1,1]=t;
                qflow_results[t+1,2*br]=MVAbase*(-modbr[br]*sin(phibr[br])*(cc[n,n,t]-cc[n,m,t])+modbr[br]*cos(phibr[br])*ss[n,m,t]);
                qflow_results[t+1,2*br+1]=MVAbase*(-modbr[br]*sin(phibr[br])*(cc[m,m,t]-cc[m,n,t])+modbr[br]*cos(phibr[br])*ss[m,n,t]);
                qlosses[br,t]=qflow_results[t+1,2*br]+qflow_results[t+1,2*br+1];
            end
        end
        qflow_results_df = DataFrame(round.(qflow_results[2:end,:], digits=digits_round), :auto);
        rename!(qflow_results_df, Symbol.(qflow_results[1,:]));
        CSV.write(string(opffolder,"/","power flow branches(MVAR).csv"), qflow_results_df);

        #Difference of angles
        anglediff_results=Array{Any}(nothing,nper+1,nbranches+1);
        anglediff_results[1,1]=string("periods");
        for br in BR, t in NPER
            n=mapbranch[br,2];
            m=mapbranch[br,3];    
            anglediff_results[1,br+1]=string("id",branch[br,:branch_id]," ",n,"-",m);
            anglediff_results[t+1,1]=t;
            anglediff_results[t+1,br+1]=-180/3.1416*atan(ss[n,m,t]/(cc[n,m,t]));
        end
        anglediff_results_df = DataFrame(round.(anglediff_results[2:end,:], digits=digits_round), :auto);
        rename!(anglediff_results_df, Symbol.(anglediff_results[1,:]));
        CSV.write(string(opffolder,"/","angle difference(degree).csv"), anglediff_results_df);

        #Active power losses across a transmission line
        plosses_results=Array{Any}(nothing,nper+1,nbranches+1);
        plosses_results[1,1]=string("periods");
        for br in BR, t in NPER
            plosses_results[1,br+1]=string("id",branch[br,:branch_id]," ",mapbranch[br,2],"-",mapbranch[br,3]);
            plosses_results[t+1,1]=t;
            plosses_results[t+1,br+1]=plosses[br,t];
        end
        plosses_results_df = DataFrame(round.(plosses_results[2:end,:], digits=digits_round), :auto);
        rename!(plosses_results_df, Symbol.(plosses_results[1,:]));
        CSV.write(string(opffolder,"/","losses(MW).csv"), plosses_results_df);

        #Reactive power losses across a transmission line
        qlosses_results=Array{Any}(nothing,nper+1,nbranches+1);
        qlosses_results[1,1]=string("periods");
        for br in BR, t in NPER
            qlosses_results[1,br+1]=string("id",branch[br,:branch_id]," ",mapbranch[br,2],"-",mapbranch[br,3]);
            qlosses_results[t+1,1]=t;
            qlosses_results[t+1,br+1]=qlosses[br,t];
        end
        qlosses_results_df = DataFrame(round.(qlosses_results[2:end,:], digits=digits_round), :auto);
        rename!(qlosses_results_df, Symbol.(qlosses_results[1,:]));
        CSV.write(string(opffolder,"/","losses(MVAR).csv"), qlosses_results_df);

        #Verifying binding constraint (only for SOCP)
        cbinding_results=Array{Any}(nothing,nper+1,nbranches+1);
        cbinding_results[1,1]=string("periods");
        for br in BR, t in NPER
            n=mapbranch[br,2];
            m=mapbranch[br,3];    
            cbinding_results[1,br+1]=string("id",branch[br,:branch_id]," ",mapbranch[br,2],"-",mapbranch[br,3]);
            cbinding_results[t+1,1]=t;
            cbinding_results[t+1,br+1]=cc[n,n,t]*cc[m,m,t]-cc[n,m,t]^2-ss[n,m,t]^2;
        end
        cbinding_results_df = DataFrame(round.(cbinding_results[2:end,:], digits=digits_round), :auto);
        rename!(cbinding_results_df, Symbol.(cbinding_results[1,:]));
        CSV.write(string(opffolder,"/","conic_constraint.csv"), cbinding_results_df);
      
        #Results of active power balance
        activbal_results=Array{Any}(nothing,nper+1,9);
        activbal_results[1,1]="Periods";
        activbal_results[1,2]="Thermal generation(MW)";
        activbal_results[1,3]="Renewable generation(MW)";
        activbal_results[1,4]="BESS discharging(MW)";
        activbal_results[1,5]="Demand response(MW)";
        activbal_results[1,6]="Load(MW)";
        activbal_results[1,7]="BESS charging(MW)";
        activbal_results[1,8]="Losses(MW)";
        activbal_results[1,9]="Error";
        for t in 1:nper 
            activbal_results[t+1,1]=t;
            activbal_results[t+1,2]=sum(pthh[g,t] for g in TH);
            activbal_results[t+1,3]=ngennodis!=0 ? sum(pree[g,t] for g in RE) : 0;
            activbal_results[t+1,4]=sum(pdiss[b,t] for b in BESS);
            activbal_results[t+1,5]=sum(pshedd[n,t] for n in N);
            activbal_results[t+1,6]=sum(pcc[n,t] for n in N);
            activbal_results[t+1,7]=sum(pchh[b,t] for b in BESS);
            activbal_results[t+1,8]=sum(plosses[br,t] for br in BR);
            activbal_results[t+1,9]=activbal_results[t+1,2]+activbal_results[t+1,3]+activbal_results[t+1,4]+activbal_results[t+1,5]+
                                    -activbal_results[t+1,6]-activbal_results[t+1,7]-activbal_results[t+1,8];
        end
        activbal_results_df = DataFrame(round.(activbal_results[2:end,:], digits=digits_round), :auto);
        rename!(activbal_results_df, Symbol.(activbal_results[1,:]));
        CSV.write(string(opffolder,"/","active power balance(MW).csv"), activbal_results_df);

        #Results of active power balance per bus
        activbalb_results=Array{Any}(nothing,nbuses*nper+1,10);
        activbalb_results[1,1]="Bus";
        activbalb_results[1,2]="Periods";
        activbalb_results[1,3]="Thermal generation(MW)";
        activbalb_results[1,4]="Renewable generation(MW)";
        activbalb_results[1,5]="BESS discharging(MW)";
        activbalb_results[1,6]="Demand response(MW)";
        activbalb_results[1,7]="Load(MW)";
        activbalb_results[1,8]="BESS charging(MW)";
        activbalb_results[1,9]="Power flow delivered(MW)";
        activbalb_results[1,10]="Error";
        for n in N
          for t in NPER
            activbalb_results[(n-1)*nper+t+1,1]=n;    
            activbalb_results[(n-1)*nper+t+1,2]=t;
            activbalb_results[(n-1)*nper+t+1,3]=sum( dis_units[g, :bus_agg]==bus[n,:bus_agg] ? pthh[g,t] : 0 for g in TH);
            activbalb_results[(n-1)*nper+t+1,4]=ngennodis!=0 ? sum( nondis_units[g, :bus_agg]==bus[n,:bus_agg] ? pree[g,t] : 0 for g in RE) : 0;
            activbalb_results[(n-1)*nper+t+1,5]=sum( bess[b, :bus_agg]==bus[n,:bus_agg] ? pdiss[b,t] : 0 for b in BESS) ;
            activbalb_results[(n-1)*nper+t+1,6]=pshedd[n,t];
            activbalb_results[(n-1)*nper+t+1,7]=pcc[n,t];
            activbalb_results[(n-1)*nper+t+1,8]=sum( bess[b, :bus_agg]==bus[n,:bus_agg] ? pchh[b,t] : 0 for b in BESS) ;
            activbalb_results[(n-1)*nper+t+1,9]=MVAbase*cc[n,n,t]*G[n,n] + MVAbase*sum( n!=m ? (cc[n,m,t]*G[n,m] - ss[n,m,t]*B[n,m]) : 0 for m in N)
            activbalb_results[(n-1)*nper+t+1,10]=activbalb_results[(n-1)*nper+t+1,3]
                                                +activbalb_results[(n-1)*nper+t+1,4]
                                                +activbalb_results[(n-1)*nper+t+1,5]
                                                +activbalb_results[(n-1)*nper+t+1,6]
                                                -activbalb_results[(n-1)*nper+t+1,7]
                                                -activbalb_results[(n-1)*nper+t+1,8]
                                                -activbalb_results[(n-1)*nper+t+1,9];
          end
        end
        activbalb_results_df = DataFrame(round.(activbalb_results[2:end,:], digits=digits_round), :auto);
        rename!(activbalb_results_df, Symbol.(activbalb_results[1,:]));
        CSV.write(string(opffolder,"/","active power balance per bus(MW).csv"), activbalb_results_df);

        #Results of reactive power balance
        reactbal_results=Array{Any}(nothing,nper+1,9);
        reactbal_results[1,1]="Periods";
        reactbal_results[1,2]="Thermal generation(MVAR)";
        reactbal_results[1,3]="Renewable generation(MVAR)";
        reactbal_results[1,4]="BESS discharging(MVAR)";
        reactbal_results[1,5]="Demand response(MVAR)";
        reactbal_results[1,6]="Load(MVAR)";
        reactbal_results[1,7]="BESS charging(MVAR)";
        reactbal_results[1,8]="Losses(MVAR)";
        reactbal_results[1,9]="Error";
        for t in 1:nper 
            reactbal_results[t+1,1]=t;
            reactbal_results[t+1,2]=sum(qthh[g,t] for g in TH);
            reactbal_results[t+1,3]=0;
            reactbal_results[t+1,4]=0;
            reactbal_results[t+1,5]=sum(qshedd[n,t] for n in N);
            reactbal_results[t+1,6]=sum(qcc[n,t] for n in N);
            reactbal_results[t+1,7]=0;
            reactbal_results[t+1,8]=sum(qlosses[br,t] for br in BR);
            reactbal_results[t+1,9]=reactbal_results[t+1,2]+reactbal_results[t+1,3]+reactbal_results[t+1,4]+reactbal_results[t+1,5]+
                                    -reactbal_results[t+1,6]-reactbal_results[t+1,7]-reactbal_results[t+1,8];
        end
        reactbal_results_df = DataFrame(round.(reactbal_results[2:end,:], digits=digits_round), :auto);
        rename!(reactbal_results_df, Symbol.(reactbal_results[1,:]));
        CSV.write(string(opffolder,"/","reactive power balance(MVAR).csv"), reactbal_results_df);

        #Results of reactive power balance per bus
        reactivbalb_results=Array{Any}(nothing,nbuses*nper+1,10);
        reactivbalb_results[1,1]="Bus";
        reactivbalb_results[1,2]="Periods";
        reactivbalb_results[1,3]="Thermal generation(MVAR)";
        reactivbalb_results[1,4]="Renewable generation(MVAR)";
        reactivbalb_results[1,5]="BESS discharging(MVAR)";
        reactivbalb_results[1,6]="Demand response(MVAR)";
        reactivbalb_results[1,7]="Load(MVAR)";
        reactivbalb_results[1,8]="BESS charging(MVAR)";
        reactivbalb_results[1,9]="Power flow delivered(MVAR)";
        reactivbalb_results[1,10]="Error";
        for n in N
          for t in NPER
            reactivbalb_results[(n-1)*nper+t+1,1]=n;    
            reactivbalb_results[(n-1)*nper+t+1,2]=t;
            reactivbalb_results[(n-1)*nper+t+1,3]=sum( dis_units[g, :bus_agg]==bus[n,:bus_agg] ? qthh[g,t] : 0 for g in TH);
            reactivbalb_results[(n-1)*nper+t+1,4]=0;
            reactivbalb_results[(n-1)*nper+t+1,5]=0;
            reactivbalb_results[(n-1)*nper+t+1,6]=qshedd[n,t];
            reactivbalb_results[(n-1)*nper+t+1,7]=qcc[n,t];
            reactivbalb_results[(n-1)*nper+t+1,8]=0;
            reactivbalb_results[(n-1)*nper+t+1,9]=-MVAbase*cc[n,n,t]*B[n,n] - MVAbase*sum( n!=m ? (cc[n,m,t]*B[n,m] + ss[n,m,t]*G[n,m]) : 0 for m in N);
            reactivbalb_results[(n-1)*nper+t+1,10]=reactivbalb_results[(n-1)*nper+t+1,3]
                                                  +reactivbalb_results[(n-1)*nper+t+1,4]
                                                  +reactivbalb_results[(n-1)*nper+t+1,5]
                                                  +reactivbalb_results[(n-1)*nper+t+1,6]
                                                  -reactivbalb_results[(n-1)*nper+t+1,7]
                                                  -reactivbalb_results[(n-1)*nper+t+1,8]
                                                  -reactivbalb_results[(n-1)*nper+t+1,9];
          end
        end
        reactivbalb_results_df = DataFrame(round.(reactivbalb_results[2:end,:], digits=digits_round), :auto);
        rename!(reactivbalb_results_df, Symbol.(reactivbalb_results[1,:]));
        CSV.write(string(opffolder,"/","reactive power balance per bus(MVAR).csv"), reactivbalb_results_df);

        #Results of LMP (Active power)
        lmpp_results=Array{Any}(nothing,nper+1,nbuses+1);
        lmpp_results[1,1]=string("periods");
        for n in N, t in NPER
                lmpp_results[1,n+1]=string("bus ", n);
                lmpp_results[t+1,[1 n+1]]=[t lmpp[n,t]];
        end
        lmpp_results_df = DataFrame(round.(lmpp_results[2:end,:], digits=digits_round), :auto);
        rename!(lmpp_results_df, Symbol.(lmpp_results[1,:]));
        CSV.write(string(opffolder,"/","LMP(Dolar per MWh).csv"), lmpp_results_df);

        #Results of LMP (Reactive power)
        lmpq_results=Array{Any}(nothing,nper+1,nbuses+1);
        lmpq_results[1,1]=string("periods");
        for n in N, t in NPER
                lmpq_results[1,n+1]=string("bus ", n);
                lmpq_results[t+1,[1 n+1]]=[t lmpq[n,t]];
        end
        lmpq_results_df = DataFrame(round.(lmpq_results[2:end,:], digits=digits_round), :auto);
        rename!(lmpq_results_df, Symbol.(lmpq_results[1,:]));
        CSV.write(string(opffolder,"/","LMP(Dolar per MVARh).csv"), lmpq_results_df);

        if consider_ctrans==1
        #Results of dual of current constraint 
        dualcurrent2_results=Array{Any}(nothing,nper+1,nbranches+1);
        dualcurrent2_results[1,1]=string("periods");
        for br in BR, t in NPER
            n=mapbranch[br,2];
            m=mapbranch[br,3];    
            dualcurrent2_results[1,br+1]=string("id",branch[br,:branch_id]," ",mapbranch[br,2],"-",mapbranch[br,3]);
            dualcurrent2_results[t+1,1]=t;
            dualcurrent2_results[t+1,br+1]=dualcurr[br,t];
        end
        dualcurrent2_results_df = DataFrame(round.( dualcurrent2_results[2:end,:], digits=digits_round), :auto);
        rename!(dualcurrent2_results_df, Symbol.(dualcurrent2_results[1,:]));
        CSV.write(string(opffolder,"/","dual_current2.csv"), dualcurrent2_results_df);
        end
 
        #Results of Objective Function
        OF_results=Array{Any}(nothing,11,2);
        OF_results[1,1]="Parameter";
        OF_results[1,2]="Value";
        OF_results[2,1]="Power shedding cost (USD)";
        OF_results[2,2]=price_p_shed*deltaT*sum(pshedd[n,t] for n in N, t in NPER);
        OF_results[3,1]="Thermal variable cost (USD)";
        OF_results[3,2]=deltaT*sum(gencost(g)*pthh[g,t] for g in TH,t in NPER);   
        OF_results[4,1]="Discharging power cost (USD)";
        OF_results[4,2]=deltaT*sum(bess[b, :DischCost]*pdiss[b,t] for b in BESS,t in NPER);
        OF_results[5,1]="Charging power cost (USD)";
        OF_results[5,2]=(-1)*deltaT*sum(bess[b, :ChargCost]*pchh[b,t] for b in BESS,t in NPER);
        OF_results[6,1]="Total cost (USD)";
        OF_results[6,2]=OF_results[2,2]+OF_results[3,2]+OF_results[4,2]+OF_results[5,2];
        OF_results[7,1]="Status";
        OF_results[7,2]=string(status1," ",status2);
        OF_results[8,1]="Type of BESS modeling";
        OF_results[8,2]=bess_modeling==1 ? "Accurate" : "Traditional";
        OF_results[9,1]="Capacity factors and load data from";
        OF_results[9,2]=string(daterun);
        OF_results[10,1]="Case name";
        OF_results[10,2]=string(case_name);
        OF_results[11,1]="Formulation";
        OF_results[11,2]=string(formulation);        
        OF_results_df = DataFrame(OF_results[2:end,:], :auto);
        rename!(OF_results_df, Symbol.(OF_results[1,:]));
        CSV.write(string(opffolder,"/","objective function.csv"), OF_results_df);

        resource=unique(plant[:,:type]);
        if consider_bess==1
        resource=vcat(resource,"_discharge_bess");
        resource=vcat(resource,"_charge_bess");
        end
        resource=vcat(resource, "curtailment");

        nresource=size(resource,1);

        genres_results=Array{Any}(nothing,nper*nresource+1,3);
        genres_results[1,1]=string("resource");
        genres_results[1,2]=string("hour");
        genres_results[1,3]=string("generation");
        for t in NPER
                for r in 1:nresource
                        genres_results[(t-1)*nresource+r+1,1]=resource[r];
                        genres_results[(t-1)*nresource+r+1,2]=t;
                        genres_results[(t-1)*nresource+r+1,3]=sum(dis_units[g,:type]==resource[r] ? pthh[g,t] : 0 for g in TH)+
                                                              (ngennodis!=0 ? sum(nondis_units[g,:type]==resource[r] ? pree[g,t] : 0 for g in RE) : 0) +
                                                              sum(resource[r]=="_discharge_bess" ? pdiss[b,t] : 0 for b in BESS)+
                                                              sum(resource[r]=="_charge_bess" ? pchh[b,t] : 0 for b in BESS)   +
                                                              (ngennodis!=0 ? sum(resource[r]=="curtailment" ? (premax[g,t]-pree[g,t]) : 0 for g in RE) : 0 );
                end
        end
        genres_results[2:end,3]=round.(genres_results[2:end,3], digits=digits_round);
        genres_results_df = DataFrame(genres_results[2:end,:], :auto);
        rename!(genres_results_df, Symbol.(genres_results[1,:]));
        CSV.write(string(opffolder,"/","generation by resource(MW).csv"), genres_results_df);

        genres_results_df[!,:resource] = convert.(String, genres_results_df[!,:resource]);
        genres_results_df[!,:hour] = convert.(Int, genres_results_df[!,:hour]);
        genres_results_df[!,:generation] = convert.(Float64,  genres_results_df[!,:generation]);        

        #Results OF in details
        OFd_results=Array{Any}(nothing,nper+1,7);
        OFd_results[1,1]=string("Periods");
        OFd_results[1,2]=string("Thermal variable cost");
        OFd_results[1,3]=string("Discharging power cost");
        OFd_results[1,4]=string("Charging power cost");
        OFd_results[1,5]=string("P demand response");
        OFd_results[1,6]=string("Q demand response");
        OFd_results[1,7]=string("Total Cost");

        for t in NPER
                OFd_results[t+1,1]=t;
                OFd_results[t+1,2]=deltaT*sum(gencost(g)*pthh[g,t] for g in TH);
                OFd_results[t+1,3]=deltaT*sum(bess[b, :DischCost]*pdiss[b,t] for b in BESS);
                OFd_results[t+1,4]=(-1)*deltaT*sum(bess[b, :ChargCost]*pchh[b,t] for b in BESS);
                OFd_results[t+1,5]=price_p_shed*deltaT*sum(pshedd[n,t] for n in N if bus[n, :Pd]!=0);
                OFd_results[t+1,6]=price_q_shed*deltaT*sum(qshedd[n,t] for n in N if bus[n, :Pd]!=0);
                OFd_results[t+1,7]=OFd_results[t+1,2]+OFd_results[t+1,3]+OFd_results[t+1,4]+OFd_results[t+1,5]+OFd_results[t+1,6];   
        end
    
        OFd_results_df = DataFrame(round.(OFd_results[2:end,:], digits=digits_round), :auto);
        rename!(OFd_results_df, Symbol.(OFd_results[1,:]));
        CSV.write(string(opffolder,"/","hourly objective function.csv"), OFd_results_df);
  
        if consider_bess==1
        #Results of bess discharging power
        dis_results=Array{Any}(nothing,nper+1,nbess+1);
        dis_results[1,1]=string("periods");
        for b in BESS, t in NPER
                 dis_results[1,b+1]=string("bess ", bess[b, :bess_id]);
                 dis_results[t+1,[1 b+1]]=[t pdiss[b,t]];
        end
        dis_results_df = DataFrame(round.(dis_results[2:end,:], digits=digits_round), :auto);
        rename!( dis_results_df, Symbol.(dis_results[1,:]));
        CSV.write(string(opffolder,"/","bess discharging power(MW).csv"), dis_results_df);

        #Results of bess charging power
        cha_results=Array{Any}(nothing,nper+1,nbess+1);
        cha_results[1,1]=string("periods");
        for b in BESS, t in NPER
                cha_results[1,b+1]=string("bess ", bess[b, :bess_id]);
                cha_results[t+1,[1 b+1]]=[t pchh[b,t]];
        end
        cha_results_df = DataFrame(round.(cha_results[2:end,:],digits=digits_round), :auto);
        rename!( cha_results_df, Symbol.(cha_results[1,:]));
        CSV.write(string(opffolder,"/","bess charging power(MW).csv"), cha_results_df);

        #Verifying binding constraint of the battery (only for SOCP)
        cbindingbess_results=Array{Any}(nothing,nper+1,nbess+1);
        cbindingbess_results[1,1]=string("periods");
        for b in BESS, t in NPER
            n=bess[b,:bus_agg][1];
            cbindingbess_results[1,b+1]=string("bess-bus ", bess[b, :bess_id], " - ", n);
            cbindingbess_results[t+1,1]=t;
            #cbindingbess_results[t+1,b+1]=(cc[n,n,t]-ibb[b,t])^2+(2*pchh[b,t]/MVAbase)^2+(2*pdiss[b,t]/MVAbase)^2-(cc[n,n,t] +  ibb[b,t])^2;
            cbindingbess_results[t+1,b+1]=(pchh[b,t]/MVAbase+pdiss[b,t]/MVAbase)^2/cc[n,n,t]-ibb[b,t];
            #cbindingbess_results[t+1,b+1]=(cc[n,n,t]*ibb[b,t]-(pchh[b,t]/MVAbase+pdiss[b,t]/MVAbase)^2)/(cc[n,n,t]*ibb[b,t]);
        end
        cbindingbess_results_df = DataFrame(round.(cbindingbess_results[2:end,:], digits=digits_round), :auto);
        rename!(cbindingbess_results_df, Symbol.(cbindingbess_results[1,:]));
        CSV.write(string(opffolder,"/","conicbess_constraint.csv"), cbindingbess_results_df);

        #Results of bess charging status
#=
        ucc_results=Array{Any}(nothing,nper+1,nbess+1);
        ucc_results[1,1]=string("periods");
        for b in BESS, t in NPER
                ucc_results[1,b+1]=string("bess ", bess[b, :bess_id]);
                ucc_results[t+1,[1 b+1]]=[t (ucc[b,t]>=0.99 ? 1 : 0)];
        end
        ucc_results_df = DataFrame(ucc_results[2:end,:], :auto);
        rename!( ucc_results_df, Symbol.(ucc_results[1,:]));
        CSV.write(string(opffolder,"/","bess charging status.csv"), ucc_results_df);

        #Results of bess discharging status
        udd_results=Array{Any}(nothing,nper+1,nbess+1);
        udd_results[1,1]=string("periods");
        for b in BESS, t in NPER
                udd_results[1,b+1]=string("bess ", bess[b, :bess_id]);
                udd_results[t+1,[1 b+1]]=[t (udd[b,t]>=0.99 ? 1 : 0)]   ;
        end
        udd_results_df = DataFrame(udd_results[2:end,:], :auto);
        rename!( udd_results_df, Symbol.(udd_results[1,:]));
        CSV.write(string(opffolder,"/","bess discharging status.csv"), udd_results_df);
=#        
        #Results of bess energy
        e_results=Array{Any}(nothing,nper+1,nbess+1);
        e_results[1,1]=string("periods");
        for b in BESS, t in NPER
                e_results[1,b+1]=string("bess ", bess[b, :bess_id]);
                e_results[t+1,[1 b+1]]=[t ee[b,t]];
        end
        e_results_df = DataFrame(round.(e_results[2:end,:],digits=digits_round), :auto);
        rename!( e_results_df, Symbol.(e_results[1,:]));
        CSV.write(string(opffolder,"/","bess energy (MWh).csv"), e_results_df);
        end

        println(string("All CSV files were printed successfully"))