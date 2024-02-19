 #File for printing graphs
       
        if consider_bess==1
        sample_battery=1;
        sample_bess = bar(1:1:24, 
                    [pdiss.data[sample_battery,:] pchh.data[sample_battery,:] ],
                    label=["Discharging power" "Charging power"], 
                #    title=string("Operation of the battery in bus ", bess[sample_battery, :bus_agg]), 
                    xticks=1:1:nper,
                    xlabel="Hour (h)", 
                    ylabel="Power (MW)",
                    legend=:bottomright, 
                    linewidth=1,
                    left_margin = 0.5Plots.mm, right_margin = 12Plots.mm,
                    bottom_margin=3Plots.mm,
                    titlefontsize=11,
                    labelfontsize=11,
                    xtickfontsize=9,
                    ytickfontsize=9,
                    legendfontsize=7);

        sample_bess = plot!(twinx(),
                    ee.data[sample_battery,:],
                    color=:green,
                    xticks = :false,   
                    left_margin = 0.5Plots.mm, right_margin = 14Plots.mm,
                    bottom_margin=3Plots.mm,
                    label="SOC",
                    ylabel="Energy (MWh)",
                    grid=true,
                    grid_linewidth=1,      
                    linewidth=2,
                    legend=:topright,
                    titlefontsize=11,
                    labelfontsize=11,
                    xtickfontsize=7,
                    ytickfontsize=9,
                    legendfontsize=6);

        namegraphic = string(opffolder,"/","Operation of a battery.pdf") 
        savefig(sample_bess, namegraphic);
        end
        #Print the generation by resource

        # Rescale hours (to start from zero)
        genres_results_df.hour = genres_results_df.hour;

        gen_res=genres_results_df |>
        @vlplot(:area, 
        x={:hour, title="Hour (h)", axis={values=0:2:24}, scale={ domain =[0, 24]}}, y={:generation, title="Power (MW)", stack=:zero}, 
        color={"resource:n", scale={scheme="category10"}})
        save(string(opffolder,"/","generation by resource(MW).pdf"),  gen_res)
      
        println(string("All plots were printed successfully"))
       