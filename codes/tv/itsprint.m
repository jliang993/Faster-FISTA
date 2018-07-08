function itsprint(str, k)
%ITSPRINT displays the iteration information

if(k==1)
    fprintf(str);
end

fprintf('%s%s', char(8*ones(1,length(str))), str);