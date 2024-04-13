%% Plot component data function
% Plot ASM2d-N2O 24 components in ydata time-series, ydata can be repeated

function plot_comp(t, ydata, opt)
    arguments
        t(:,:) double {mustBeNumeric};
    end
    arguments (Repeating)
        ydata (:,24) double {mustBeNumeric};
    end
    arguments
        opt.f = figure;
        opt.numRow (1,1) {mustBeInteger, mustBePositive} = 4;
        opt.numCol (1,1) {mustBeInteger, mustBePositive} = 6;
        opt.yLabels (1,:) string = {'True data'}; % must be paired with the input ydata
        opt.yName (1,24) string = {"S_{O_2}","S_F","S_A","S_{NH_4}","S_{NH_2OH}","S_{N_2O}", ...
            "S_{NO}","S_{NO_2}","S_{NO_3}","S_{PO_4}","S_I","S_{ALK}","N_2","X_I","X_S","X_H", ...
            "X_{PAO}","X_{PP}","X_{PHA}","X_{AOB}","X_{NOB}","X_{TSS}","X_{MeOH}","X_{Mep}"};
        opt.yUnit (1,24) string = {"g(O_2)/m^3","g(COD)/m^3","g(COD)/m^3","g(N)/m^3","g(N)/m^3","g(N)/m^3",...
            "g(N)/m^3","g(N)/m^3","g(N)/m^3","g(P)/m^3","g(COD)/m^3","mol(HCO_3^-)/m^3",...
            "g(N)/m^3","g(COD)/m^3","g(COD)/m^3","g(COD)/m^3","g(COD)/m^3","g(P)/m^3",...
            "g(COD)/m^3","g(COD)/m^3","g(COD)/m^3","g(TSS)/m^3","g(TSS)/m^3","g(TSS)/m^3"};
        opt.dyName (1,24) string = {"S_{O_2}'","S_F'","S_A'","S_{NH_4}'","S_{NH_2OH}'","S_{N_2O}'", ...
            "S_{NO}'","S_{NO_2}'","S_{NO_3}'","S_{PO_4}'","S_I'","S_{ALK}'","N_2'","X_I'","X_S'","X_H'", ...
            "X_{PAO}'","X_{PP}'","X_{PHA}'","X_{AOB}'","X_{NOB}'","X_{TSS}'","X_{MeOH}'","X_{Mep}'"};
        opt.dyUnit (1,24) string = {"g(O_2)/m^3/d","g(COD)/m^3/d","g(COD)/m^3/d","g(N)/m^3/d","g(N)/m^3/d","g(N)/m^3/d",...
            "g(N)/m^3/d","g(N)/m^3/d","g(N)/m^3/d","g(P)/m^3/d","g(COD)/m^3/d","mol(HCO_3^-)/m^3/d",...
            "g(N)/m^3/d","g(COD)/m^3/d","g(COD)/m^3/d","g(COD)/m^3/d","g(COD)/m^3/d","g(P)/m^3/d",...
            "g(COD)/m^3/d","g(COD)/m^3/d","g(COD)/m^3/d","g(TSS)/m^3/d","g(TSS)/m^3/d","g(TSS)/m^3/d"};
        opt.figureTitle string ;
        opt.isDy (1,1) logical {mustBeNumericOrLogical} = false;
        opt.showUnit (1,1) logical {mustBeNumericOrLogical} = false;
        opt.showRMSE (1,1) logical {mustBeNumericOrLogical} = false;
        opt.showMAE (1,1) logical {mustBeNumericOrLogical} = false;
        opt.showR2 (1,1) logical {mustBeNumericOrLogical} = false;
    end

    % opt.f.Position(3) = 1.5*opt.f.Position(3);
    colororder(["#0072BD" "#A2142F" "#EDB120" ]);
    % linestyleorder("mixedstyles");
    linestyles = [":","-.","-"];
    linestyleorder(linestyles);

    thour = t*24; % change from day to hour for plotting
    no_ydata = numel(ydata);

    if opt.isDy
        Name = opt.dyName;
        Unit = opt.dyUnit;
    else
        Name = opt.yName;
        Unit = opt.yUnit;
    end

    if opt.showRMSE
        if no_ydata >= 2
            % calcuate RMSE of first two inputs of ydata only            
            strRMSE = sprintf('%s %.2f ', "RMSE=", rmse(ydata{1}, ydata{2},"all")); % rmse = sqrt(mean((dyColl - dyPred).^2))
        else
            opt.showRMSE = false;
        end
    end

    if opt.showMAE
        if no_ydata >= 2
            % calcuate MAE of first two inputs of ydata only            
            strMAE = sprintf('%s %.2f ', "MAE=", mean(abs(ydata{1}-ydata{2}),"all"));
        else
            opt.showMAE = false;
        end
    end

    if opt.showR2
        if no_ydata >= 2
            % calcuate R2 of first two inputs of ydata only  
            r2 = corrcoef(ydata{1}, ydata{2});
            strR2 = sprintf('%s %.2f ', "R^2=", r2(2,1));
        else
            opt.showR2 = false;
        end
    end
      
    for k = 1:24
        subplot(opt.numRow, opt.numCol, k)
        for j = 1:no_ydata
            plot(thour, ydata{j}(:,k), LineWidth=1)
            hold on
        end
        hold off
        
        if opt.showRMSE 
            title(sprintf('%s %s %.2f ', Name(k), ' | RMSE=', rmse(ydata{1}(:,k), ydata{2}(:,k),"all")));
        else
            title(Name(k));
        end

        xlim([0 6])
        xticks(0:2:6)
        xticklabels({'0','2','4','6h'})
        ytickformat('%.1g')
        if opt.showUnit
            ylabel(Unit(k))
        end

        if k==2 
            legend(opt.yLabels,'Location','best')
        end
        grid on
    end

    % show legend in the bottom right of the figure
    % lgd = legend(opt.yLabels,'Orientation','horizontal');
    % lgd.Position(1) = 0.62;
    % lgd.Position(2) = 0.03;
    
    % show index value in the bottom right as annotation
    % if  opt.showRMSE
    %     delete(findall(gcf,'type','annotation'));
    %     anno = annotation(opt.f, 'textbox',[.1 .03 .2 .03],'String',strRMSE,'FitBoxToText','on');
    %     anno.LineStyle = "none";
    % end

    if isfield(opt,'figureTitle')
        strTitle = opt.figureTitle;
        if opt.showRMSE
            strTitle = strTitle + " | " + strRMSE;
        end
        if opt.showMAE
            strTitle = strTitle + " | " + strMAE;
        end
        if opt.showR2
            strTitle = strTitle + " | " + strR2;
        end
        sgtitle(strTitle);
    else
        if  opt.showRMSE
            strTitle = strRMSE;
            if opt.showMAE
                strTitle = strTitle + " | " + strMAE;
            end
            if opt.showR2
                strTitle = strTitle + " | " + strR2;
            end
            sgtitle(strTitle);
        else
            if opt.showMAE
                strTitle = strMAE;
                if opt.showR2
                    strTitle = strTitle + " | " + strR2;
                end
                sgtitle(strTitle);                    
            else
                if opt.showR2
                    strTitle = strR2;
                    sgtitle(strTitle); 
                end                      
            end                
        end
    end
end