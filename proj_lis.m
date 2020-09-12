function output=proj_lis(input,bit)
if bit==0
    output=input;
    return
end
if bit==-1
    output=ones(size(input));
    return
end
output=nan(size(input));
step=2*pi/(2^bit);
list=exp((-pi:step:(pi-step))*1i);
for i=1:length(input)
    [~,I]=min( abs(list-input(i)));
    output(i)=list(I);
end
end