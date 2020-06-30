clear all
close all

%Mass balances for abiotic static and transient cases
%Order of compartments 
    % For breakthrough curve estimates: Effluent, Retained
    % For direct measurement estimates: Effluent, Headspace Before,
    % Headspace After

static1btn2o = [0.9689, 0.0331, 0];
static1dmn2o = [0.9689, 0.0652, 0.0874];

static1btch4 = [0.8948, 0.1052, 0];
static1dmch4 = [0.8948, 0.0465, 0.0725];

static2btn2o = [0.9808, 0.0192, 0];
static2dmn2o = [0.9808, 0.0836, 0.0971];

static2btch4 = [0.8593, 0.1407, 0];
static2dmch4 = [0.8593, 0.0587, 0.0793];

tran1btn2o = [0.9243, 0.0757, 0];
tran1dmn2o = [0.9243, 0.0409, 0.0325];

tran1btch4 = [0.8140, 0.1860, 0];
tran1dmch4 = [0.8140, 0.0146, 0.0432];

tran2btn2o = [0.9188, 0.0812, 0];
tran2dmn2o = [0.9188, 0.0286, 0.1017];

tran2btch4 = [0.7952, 0.2048, 0];
tran2dmch4 = [0.7952, 0.0161, 0.1378];

null = [0, 0, 0]
btbarinfon2o = vertcat(static1btn2o, static1dmn2o, static2btn2o, static2dmn2o, tran1btn2o, tran1dmn2o, tran2btn2o, tran2dmn2o);
figure(2000)
% x = categorical(["Static Case 1 Breakthrough Curve MB" "Static Case 1 Direct Measurement MB" "Static Case 2 Breakthrough Curve MB" "Static Case 2 Direct Measurement MB" "Transient Case 1 Breakthrough Curve MB" "Transient Case 1 Direct Measurement MB" "Transient Case 2 Breakthrough Curve MB" "Transient Case 2 Direct Measurement MB"]);
% x = reordercats(x,{'Static Case 1 Breakthrough Curve MB' 'Static Case 1 Direct Measurement MB' 'Static Case 2 Breakthrough Curve MB' 'Static Case 2 Direct Measurement MB' 'Transient Case 1 Breakthrough Curve MB' 'Transient Case 1 Direct Measurement MB' 'Transient Case 2 Breakthrough Curve MB' 'Transient Case 2 Direct Measurement MB' });
xloc = [1, 1.75, 2.75, 3.5, 4.5, 5.25, 6.25, 7];
bbt = bar(xloc, btbarinfon2o, 'stacked');
hold on 
%Recolor breakthrough curve based retained mass
btcolor = [156/255 172/255 165/255]; 
effcolor = [20/255 135/255 240/255];
hsbeforecolor = 'r'
hsaftercolor = 'y'
r1 = rectangle('Position', [0.7, static1btn2o(1), 0.6, static1btn2o(2)], 'FaceColor',btcolor,'EdgeColor','k')
r2 = rectangle('Position', [2.45, static2btn2o(1), 0.6, static2btn2o(2)], 'FaceColor',btcolor,'EdgeColor','k')
r3 = rectangle('Position', [4.2, tran1btn2o(1), 0.6, tran1btn2o(2)], 'FaceColor',btcolor,'EdgeColor','k')
r4 = rectangle('Position', [5.95, tran2btn2o(1), 0.6, tran2btn2o(2)], 'FaceColor',btcolor,'EdgeColor','k')
bbt(1).FaceColor = effcolor;
bbt(2).FaceColor = hsbeforecolor;
bbt(3).FaceColor = hsaftercolor;
%legend('N_2O', 'SF_6', 'C_2H_6', 'CH_4', 'He', 'Location', 'northwest')
ylabel('Normalized Mass [M/M_T]')
set(gca,'xtick',[])
xticks(xloc)
xticklabels({'BC\newline \newline p' 'DM', 'BC', 'DM', 'BC', 'DM', 'BC', 'DM'})
%xtickangle(45)
xlim([0.5 7.5])
set(gca,'TickDir','out');
set(gca, 'FontSize', 22);
text(1.375,-0.16,'Static Case 1','FontSize', 22, 'HorizontalAlignment', 'center')
text(3.125,-0.16,'Static Case 2','FontSize', 22, 'HorizontalAlignment', 'center')
text(4.875,-0.16,'Transient Case 1','FontSize', 22, 'HorizontalAlignment', 'center')
text(6.625,-0.16,'Transient Case 2','FontSize', 22, 'HorizontalAlignment', 'center')


btbarinfoch4 = vertcat(static1btch4, static1dmch4, static2btch4, static2dmch4, tran1btch4, tran1dmch4, tran2btch4, tran2dmch4);
figure(4000)
% x = categorical(["Static Case 1 Breakthrough Curve MB" "Static Case 1 Direct Measurement MB" "Static Case 2 Breakthrough Curve MB" "Static Case 2 Direct Measurement MB" "Transient Case 1 Breakthrough Curve MB" "Transient Case 1 Direct Measurement MB" "Transient Case 2 Breakthrough Curve MB" "Transient Case 2 Direct Measurement MB"]);
% x = reordercats(x,{'Static Case 1 Breakthrough Curve MB' 'Static Case 1 Direct Measurement MB' 'Static Case 2 Breakthrough Curve MB' 'Static Case 2 Direct Measurement MB' 'Transient Case 1 Breakthrough Curve MB' 'Transient Case 1 Direct Measurement MB' 'Transient Case 2 Breakthrough Curve MB' 'Transient Case 2 Direct Measurement MB' });
xloc = [1, 1.75, 2.75, 3.5, 4.5, 5.25, 6.25, 7];
bbt = bar(xloc, btbarinfoch4, 'stacked');
hold on 
%Recolor breakthrough curve based retained mass
btcolor = [156/255 172/255 165/255]; 
effcolor = [20/255 135/255 240/255];
hsbeforecolor = 'r'
hsaftercolor = 'y'
r1 = rectangle('Position', [0.7, static1btch4(1), 0.6, static1btch4(2)], 'FaceColor',btcolor,'EdgeColor','k')
r2 = rectangle('Position', [2.45, static2btch4(1), 0.6, static2btch4(2)], 'FaceColor',btcolor,'EdgeColor','k')
r3 = rectangle('Position', [4.2, tran1btch4(1), 0.6, tran1btch4(2)], 'FaceColor',btcolor,'EdgeColor','k')
r4 = rectangle('Position', [5.95, tran2btch4(1), 0.6, tran2btch4(2)], 'FaceColor',btcolor,'EdgeColor','k')
bbt(1).FaceColor = effcolor;
bbt(2).FaceColor = hsbeforecolor;
bbt(3).FaceColor = hsaftercolor;
ylabel('Normalized Mass [M/M_T]')
set(gca,'xtick',[])
xticks(xloc)
xticklabels({'BC\newline \newline p' 'DM', 'BC', 'DM', 'BC', 'DM', 'BC', 'DM'});
%xtickangle(45)
xlim([0.5 7.5])
set(gca,'TickDir','out');
set(gca, 'FontSize', 22);
text(1.375,-0.16,'Static Case 1','FontSize', 22, 'HorizontalAlignment', 'center')
text(3.125,-0.16,'Static Case 2','FontSize', 22, 'HorizontalAlignment', 'center')
text(4.875,-0.16,'Transient Case 1','FontSize', 22, 'HorizontalAlignment', 'center')
text(6.625,-0.16,'Transient Case 2','FontSize', 22, 'HorizontalAlignment', 'center')

%I am legend
figure(2007)
legenddata1 = [1, 1, 1, 1];
legenddata2 = [2, 2, 2, 2];
legenddatacomb = vertcat(legenddata1, legenddata2);
legendbar = bar (legenddatacomb, 'stacked');
legendbar(1).FaceColor = effcolor;
legendbar(2).FaceColor = btcolor;
legendbar(3).FaceColor = hsbeforecolor;
legendbar(4).FaceColor = hsaftercolor;
lgd = legend('Effluent Mass__', 'Modeled Retained Mass__', 'Measured HS Evolution__', 'Measured Post-Drainage Bubble Release__', 'Location', 'northwest', 'Orientation','horizontal')
lgd.FontSize = 20;
% pos=lgd.Position;        
% pos(3)=0.75;
% pos(1)=0.25;
% lgd.Position=pos;       
ylim([0 20])

