files = dir('/archival/projects/SPLaT/data/beh/**/sub-*_task-naturalistic_beh_*.csv');
animals = {'dog','cat','frog','spider'}; loom = {'none' 'loom'};


for f = 1:length(files);
    t = readtable([files(f).folder filesep files(f).name]);
    for a=1:length(animals)
        for l=0:1
            wh_rows = strcmp(t.animal_type,animals{a}) & t.has_loom==l;
            X(f,a,l+1,1) = mean(t.fear_rating(wh_rows));
            X(f,a,l+1,2) = mean(t.pleasantness_rating(wh_rows));
            X(f,a,l+1,3) =mean(t.arousal_rating(wh_rows));

        end
    end
end

%%
measures = {'Fear' 'Valence' 'Arousal'};
for m=1:3
    for f = 1:size(X,1)

        rX(f,:)=[squeeze(X(f,:,1,m)) squeeze(X(f,:,2,m))];
    end

    anovaTable=array2table(rX,'VariableNames',{'x1','x2','x3','x4','x5','x6','x7','x8'});

    factors = categorical([1 1 1 1 2 2 2 2 ;
        [1:4 1:4]]');
    rm = fitrm(anovaTable,"x1-x8 ~ 1",'WithinDesign',factors,'WithinModel',"orthogonalcontrasts");
    ranova(rm,'WithinModel','w1*w2')

    t= multcompare(rm,'w1','by','w2');
    t = t(1:2:end,:)

    figure(1)
    subplot(1,6,2*(m-1)+1)
    barplot_columns(rX(:,1:4),'dolines','nofig');
    set(gca,'XTickLabel',animals)
    xlabel('Looming Absent')
    ylabel(measures{m})
    ylim([-50 150])

    subplot(1,6,2*(m-1)+2)
barplot_columns(rX(:,5:end),'dolines','nofig');
    set(gca,'XTickLabel',animals)
    ylim([-50 150])
    ylabel(measures{m})
    xlabel('Looming Present')


    tbl = margmean(rm,{'w1','w2'})

    figure(1+m); 
    barplot_columns(rX(:,5:end)-rX(:,1:4),'nofig')
    ylabel(['\Delta ' measures{m}])
    xlabel('Effect of Looming')
    set(gca,'XTickLabel',animals)

end
