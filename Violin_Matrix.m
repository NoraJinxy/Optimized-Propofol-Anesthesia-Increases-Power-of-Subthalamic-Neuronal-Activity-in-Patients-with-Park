function [average,name]=Violin_Matrix(type,num1,num2,fq1,fq2,f)
if type==1
    for i=1:num1
        average(i,1)=mean(fq1(i,find(abs(f-0.5)<1e-5):find(abs(f-3)<1e-5)),2);
        name{i,1}='delta L';
        average(i,3)=mean(fq1(i,find(abs(f-4)<1e-5):find(abs(f-8)<1e-5)),2);
        name{i,3}='theta L';
        average(i,5)=mean(fq1(i,find(abs(f-8)<1e-5):find(abs(f-14)<1e-5)),2);
        name{i,5}='alpha L';
        average(i,7)=mean(fq1(i,find(abs(f-14)<1e-5):find(abs(f-31)<1e-5)),2);
        name{i,7}='beta L';
        average(i,9)=mean(fq1(i,find(abs(f-31)<1e-5):find(abs(f-100)<1e-5)),2);
        name{i,9}='gamma L';

        average(i,2)=mean(fq2(i,find(abs(f-0.5)<1e-5):find(abs(f-3)<1e-5)),2);
        name{i,2}='delta R';
        average(i,4)=mean(fq2(i,find(abs(f-4)<1e-5):find(abs(f-8)<1e-5)),2);
        name{i,4}='theta R';
        average(i,6)=mean(fq2(i,find(abs(f-8)<1e-5):find(abs(f-14)<1e-5)),2);
        name{i,6}='alpha R';
        average(i,8)=mean(fq2(i,find(abs(f-14)<1e-5):find(abs(f-31)<1e-5)),2);
        name{i,8}='beta R';
        average(i,10)=mean(fq2(i,find(abs(f-31)<1e-5):find(abs(f-100)<1e-5)),2);
        name{i,10}='gamma R';
    end
else if type==2
      for i=1:num1
        average1(i,1)=mean(fq1(i,find(abs(f-0.5)<1e-5):find(abs(f-3)<1e-5)),2);
        name1{i,1}='delta TCI';
        average1(i,2)=mean(fq1(i,find(abs(f-4)<1e-5):find(abs(f-8)<1e-5)),2);
        name1{i,2}='theta TCI';
        average1(i,3)=mean(fq1(i,find(abs(f-8)<1e-5):find(abs(f-14)<1e-5)),2);
        name1{i,3}='alpha TCI';
        average1(i,4)=mean(fq1(i,find(abs(f-14)<1e-5):find(abs(f-31)<1e-5)),2);
        name1{i,4}='beta TCI';
        average1(i,5)=mean(fq1(i,find(abs(f-31)<1e-5):find(abs(f-100)<1e-5)),2);
        name1{i,5}='gamma TCI';
    end
    for i=1:num2
        average2(i,1)=mean(fq2(i,find(abs(f-0.5)<1e-5):find(abs(f-3)<1e-5)),2);
        name2{i,1}='delta Awake';
        average2(i,2)=mean(fq2(i,find(abs(f-4)<1e-5):find(abs(f-8)<1e-5)),2);
        name2{i,2}='theta Awake';
        average2(i,3)=mean(fq2(i,find(abs(f-8)<1e-5):find(abs(f-14)<1e-5)),2);
        name2{i,3}='alpha Awake';
        average2(i,4)=mean(fq2(i,find(abs(f-14)<1e-5):find(abs(f-31)<1e-5)),2);
        name2{i,4}='beta Awake';
        average2(i,5)=mean(fq2(i,find(abs(f-31)<1e-5):find(abs(f-100)<1e-5)),2);
        name2{i,5}='gamma Awake';
    end  
    average=[average2(:);average1(:)];
    name=[name2(:);name1(:)];
    end        
end
end