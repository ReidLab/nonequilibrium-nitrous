clear all
close all

%vg/vw estimates for gases in various cases
%Order of gases - N2O, SF6, Ethane, Methane, Helium

gb1 = [0.054, 0.007, 0.028, 0.019, 0.019];
gb2 = [0.065, 0.005, 0.031, 0.018, 0.020];

ss1 = [0.107, 0.101, 0.113, 0.106, 0.140];
ss2 = [0.095, 0.096, 0.109, 0.088, 0.140];

tran1 = [0.208, 0.284, 0.231, 0.203, 0.189];
tran2 = [0.214, 0.277, 0.226, 0.200, 0.151];

bio1 = [0.336, 0.339, 0.360];%bio1 = [NaN, 0.336, 0.339, NaN, 0.360];
bio2 = [0.321, 0.369, 0.369];%bio2 = [~-1, 0.321, 0.369, ~-1, 0.369];

biojoined = [0.336,0.321; 0.339,0.369; 0.360,0.369];

barinfo = vertcat(gb1, gb2, ss1, ss2, tran1, tran2);%, bio1, bio2);

figure(1000)
x = categorical(["Glass Bead 1" "Glass Bead 2" "Static 1" "Static 2" "Transient 1" "Transient 2" "Biotic 1" "Biotic 2"]);
x = reordercats(x,{'Glass Bead 1' 'Glass Bead 2' 'Static 1' 'Static 2' 'Transient 1' 'Transient 2' 'Biotic 1' 'Biotic 2' });
b = bar(barinfo);
b(1).FaceColor = 'm';
b(2).FaceColor = 'g';
b(3).FaceColor = 'b';
b(4).FaceColor = 'c';
b(5).FaceColor = 'r';
hold on
b2_1 = bar([6.84],bio1(1), 0.125); %4.84 4.98], bio1, 0.8)
b2_2 = bar([7],bio1(2), 0.125);
b2_3 = bar([7.16],bio1(3), 0.125);
b2_1.FaceColor = 'g';
b2_2.FaceColor = 'b';
b2_3.FaceColor = 'r';
b3_1 = bar([7.84], bio2(1), 0.125);
b3_2 = bar([8], bio2(2), 0.125);
b3_3 = bar([8.16],bio2(3), 0.125);
b3_1.FaceColor = 'g';
b3_2.FaceColor = 'b';
b3_3.FaceColor = 'r';
lgd = legend('N_2O', 'SF_6', 'C_2H_6', 'CH_4', 'He', 'Location', 'northwest');
set(lgd,'FontSize',26);
ylabel('V_g/V_w [-]')
set(gca,'xtick',[])
%xlabel('Experimental Setup')
xticks([1 2 3 4 5 6 7 8])
xticklabels({' Glass\newlineBead 1' ' Glass\newlineBead 2' '    Static\newlineWoodchip 1','    Static\newlineWoodchip 2','  Transient\newlineWoodchip 1','  Transient\newlineWoodchip 2','     Biotic\newlineWoodchip 1','     Biotic\newlineWoodchip 2'})
%xtickangle(45)
%gb1pos = [0.7 0.86 1.02 1.18 1.34]
%gb1lab = text(gb1pos,gb1+0.008,num2str(gb1'),'vert','bottom','horiz','center'); 
%set(gb1lab,'Rotation',60)
xlim([0.5 8.5])
set(gca,'TickDir','out');
set(gca, 'FontSize', 30);