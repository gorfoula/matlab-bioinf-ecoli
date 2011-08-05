function [SumsHydro Start Fin]=SumPho(Seq_Name,Seq,Peptide,outFile,save_)

[PhoMatrix Pho]=Hydrophobicity(Peptide,'K');
[PolMatrix Pol]=Polarity(Peptide,'Z');

[Row,Col,Zax]=size(PhoMatrix);

if(Col<=20)
    windend=Col;
else
    windend=20;
end

MaxSum=-1000;
MaxPol=1000;
CurWindow=0;
%windend=round(0.8*Col)

for window=4:1:windend  %% 3.6 residues per turn 
    SumsHydro=conv(ones(Row,window,Zax),PhoMatrix);  %% all sums of Hydrophobicity for Window from 9 to 20
    SumsPolar=conv(ones(Row,window,Zax),PolMatrix);
        
    if(MaxPol>min(SumsPolar))
        CurSumsPolar=SumsPolar;
        MaxPol=min(SumsPolar);
        posPolMax=find((SumsPolar==MaxPol),1,'first');
        [window posPolMax MaxPol]
    end
    
    if(MaxSum<max(SumsHydro))
        CurWindow=window;
        CurSumsHydro=SumsHydro;
        MaxSum=max(SumsHydro);
        posMax=find((SumsHydro==MaxSum),1,'first');
        %[CurWindow posMax MaxSum]
    end
%     figure(4);
%     plot(SumsPolar,'-gs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize',2);
%     hold on;
%     plot(PolMatrix,'--mo','LineWidth',2);        
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     plot(SumsHydro,'-bs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize',2);
%     hold on;
%     plot(PhoMatrix,'--ro','LineWidth',2);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     xlabel('Window - Position of sum');
%     ylabel('Hydrophobicity');
%     legend('Sum Pol','Pol at each position','Sum Pho','Pho at each position');
%     title([Seq_Name,' ',int2str(window)]);
%     hold off;
%     pause;
end

if(posMax<Col)
    Fin=posMax;
else
    Fin=Col;
end
Start=Fin-CurWindow+1;
Start=Start+1;
Fin=Fin+1;
[Start Fin]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3=figure(3);
plot(CurSumsPolar,'-gs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize',2);
hold on;
plot(PolMatrix,'--mo','LineWidth',2);        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(CurSumsHydro,'-bs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize',2);
hold on;
plot(PhoMatrix,'--ro','LineWidth',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Window - Position of sum');
ylabel('Hydrophobicity');
legend('Sum Pol','Pol at each position','Sum Pho','Pho at each position');
title([Seq_Name,' ',Seq(1:Start-1),'  ',Seq(Start:Fin),'  ',Seq(Fin+1:Col)]);
hold off;
pause;

    if(save_==1)
        saveas(f3,[outFile,'_',Seq_Name,'.jpg'],'jpg');
    end 

end