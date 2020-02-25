function fraction = Compute_weight_spectrum(NetFile, Start, Stop)


[Num, Den] = mason(NetFile, Start, Stop);
Num = str2sym(Num);
Den = str2sym(Den);

disp(['Start node: ',num2str(Start),' End node: ',num2str(Stop),' Weight enumerating function: ']);
fraction = simplifyFraction(Num/Den);