%[status, result] = system('get_free_nodes gpuonly')
%parts = regexp(result,'\n','split');
%parts = parts(1:end-1)
%machines = cellstr(parts);

machines = {'cn109','cn110','cn111','cn112','cn114','cn116','cn117','cn118','cn119','cn120','cn121','cn122','cn123','cn125','cn126','cn127','cn128','cn129','cn130','cn131','cn132','cn14','cn145','cn146','cn147','cn148','cn149','cn15','cn150','cn151','cn152','cn153','cn154','cn155','cn157','cn158','cn159','cn16','cn160','cn161','cn162','cn163','cn164','cn17','cn18','cn19','cn20','cn22','cn23','cn24','cn37','cn38','cn39','cn40','cn41','cn42','cn44','cn45','cn47','cn48','cn94','cn95','cn96'}

eval( MPI_Run( 'Phyllosphere_Bootstrapno01000_3', size(machines,2), machines ))
