%[status, result] = system('get_free_nodes gpuonly')
%parts = regexp(result,'\n','split');
%parts = parts(1:end-1)
%machines = cellstr(parts);

machines = {'cn27','cn28','cn29','cn30','cn31','cn32','cn33','cn34','cn35','cn36','cn61','cn62','cn63','cn64','cn65','cn66','cn67','cn68','cn69','cn70','cn71','cn72','cn97','cn98','cn99','cn100','cn101','cn102','cn103','cn104','cn105','cn106','cn107'}


eval( MPI_Run( 'Phyllosphere_Bootstrapno01000_3', size(machines,2), machines ))
