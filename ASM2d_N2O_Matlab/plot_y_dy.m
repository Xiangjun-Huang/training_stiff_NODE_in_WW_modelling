% Plot functions
% Visualise y and y' (dy) in one figure

function plot_y_dy(t,y,dy)
    nameFeature = ["S_{O_2}","S_F","S_A","S_{NH_4}","S_{NH_2OH}","S_{N_2O}",...
                    "S_{NO}","S_{NO_2}","S_{NO_3}","S_{PO_4}","S_I","S_{ALK}",...
                    "N_2","X_I","X_S","X_H","X_{PAO}","X_{PP}",...
                    "X_{PHA}","X_{AOB}","X_{NOB}","X_{TSS}","X_{MeOH}","X_{Mep}"];
    nameUnit = ["g(O_2)/m^3","g(COD)/m^3","g(COD)/m^3","g(N)/m^3","g(N)/m^3","g(N)/m^3",...
                "g(N)/m^3","g(N)/m^3","g(N)/m^3","g(P)/m^3","g(COD)/m^3","mol(HCO_3^-)/m^3",...
                "g(N)/m^3","g(COD)/m^3","g(COD)/m^3","g(COD)/m^3","g(COD)/m^3","g(P)/m^3",...
                "g(COD)/m^3","g(COD)/m^3","g(COD)/m^3","g(TSS)/m^3","g(TSS)/m^3","g(TSS)/m^3"];  
    f=figure;
    f.Position(3) = 2*f.Position(3);
    f.Position(4) = 2.5*f.Position(4);
    for k = 1:24
        subplot(4,6,k)
        yyaxis left
        plot(t, y(:,k))
        title(nameFeature(k)) 
        xlim([0 6])
        ylim([0 inf])
        xticks(0:2:6)
        xticklabels({'0','2','4','6h'})
        ylabel(nameUnit(k))
    
        grid on
        yyaxis right
        plot(t,dy(:,k),"--")
        ylabel(nameUnit(k)+'/d')
        if k==1 
            legend('True Y','True Y''','Location','southeast')
        end
    end
    titleText = sprintf("True Y and True dY'");
    sgtitle(titleText);
end