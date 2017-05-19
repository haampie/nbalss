function bifurcation_diagram(column)

  names = {'typ', 'Re', 'Max psi', 'Min psi', 'Max + min psi', 'ps(n/8, m/8)', 'ps(3n/8, 3m/8)'};

  data_sets = {
    % *  typ   par( 5)    max psi     min psi    max + min psi    ps(n/8,m/8)    ps(3n/8,3m/8)

    % symm_re_29
    [
    0  0.163925E+02   0.109788E+01  -0.109788E+01  -0.222045E-15   0.410954E+00   0.484386E+00;
    0  0.167841E+02   0.111177E+01  -0.111177E+01   0.666134E-14   0.410923E+00   0.493089E+00;
    0  0.171746E+02   0.112572E+01  -0.112572E+01   0.621725E-14   0.410820E+00   0.502440E+00;
    0  0.175638E+02   0.113972E+01  -0.113972E+01   0.122125E-13   0.410647E+00   0.512452E+00;
    0  0.179514E+02   0.115377E+01  -0.115377E+01   0.000000E+00   0.410404E+00   0.523132E+00;
    0  0.183373E+02   0.116835E+01  -0.116835E+01   0.155431E-14   0.410094E+00   0.534484E+00;
    0  0.187213E+02   0.118321E+01  -0.118321E+01  -0.355271E-14   0.409719E+00   0.546511E+00;
    0  0.191034E+02   0.119813E+01  -0.119813E+01  -0.266454E-14   0.409280E+00   0.559213E+00;
    0  0.194832E+02   0.121310E+01  -0.121310E+01  -0.199840E-14   0.408781E+00   0.572587E+00;
    0  0.198609E+02   0.122842E+01  -0.122842E+01  -0.976996E-14   0.408223E+00   0.586626E+00;
    0  0.202361E+02   0.124424E+01  -0.124424E+01   0.732747E-14   0.407609E+00   0.601324E+00;
    0  0.206090E+02   0.126038E+01  -0.126038E+01   0.688338E-14   0.406941E+00   0.616670E+00;
    0  0.209793E+02   0.127654E+01  -0.127654E+01   0.137668E-13   0.406221E+00   0.632653E+00;
    0  0.213469E+02   0.129271E+01  -0.129271E+01  -0.113243E-13   0.405452E+00   0.649258E+00;
    0  0.217119E+02   0.130888E+01  -0.130888E+01  -0.377476E-14   0.404635E+00   0.666471E+00;
    0  0.220742E+02   0.132545E+01  -0.132545E+01  -0.244249E-14   0.403772E+00   0.684275E+00;
    0  0.224336E+02   0.134243E+01  -0.134243E+01  -0.377476E-14   0.402867E+00   0.702650E+00;
    0  0.227902E+02   0.135942E+01  -0.135942E+01  -0.199840E-14   0.401920E+00   0.721580E+00;
    0  0.231438E+02   0.137642E+01  -0.137642E+01  -0.532907E-14   0.400934E+00   0.741042E+00;
    0  0.234944E+02   0.139343E+01  -0.139343E+01   0.222045E-14   0.399911E+00   0.761017E+00;
    0  0.238420E+02   0.141071E+01  -0.141071E+01  -0.377476E-14   0.398851E+00   0.781482E+00;
    0  0.241865E+02   0.142851E+01  -0.142851E+01   0.666134E-15   0.397758E+00   0.802417E+00;
    0  0.245278E+02   0.144632E+01  -0.144632E+01  -0.421885E-14   0.396632E+00   0.823799E+00;
    0  0.248659E+02   0.146414E+01  -0.146414E+01  -0.128786E-13   0.395475E+00   0.845605E+00;
    0  0.252007E+02   0.148197E+01  -0.148197E+01  -0.131006E-13   0.394289E+00   0.867813E+00;
    0  0.255323E+02   0.149982E+01  -0.149982E+01  -0.666134E-15   0.393075E+00   0.890402E+00;
    0  0.258605E+02   0.151839E+01  -0.151839E+01   0.000000E+00   0.391835E+00   0.913348E+00;
    0  0.261853E+02   0.153700E+01  -0.153700E+01  -0.666134E-15   0.390570E+00   0.936632E+00;
    0  0.265068E+02   0.155562E+01  -0.155562E+01  -0.111022E-14   0.389282E+00   0.960231E+00;
    0  0.268248E+02   0.157426E+01  -0.157426E+01   0.266454E-14   0.387972E+00   0.984126E+00;
    0  0.271393E+02   0.159291E+01  -0.159291E+01  -0.599520E-14   0.386640E+00   0.100830E+01;
    0  0.274503E+02   0.161181E+01  -0.161181E+01   0.155431E-14   0.385290E+00   0.103272E+01;
    0  0.277577E+02   0.163120E+01  -0.163120E+01  -0.488498E-14   0.383921E+00   0.105739E+01;
    0  0.280616E+02   0.165060E+01  -0.165060E+01   0.173195E-13   0.382535E+00   0.108227E+01;
    0  0.283619E+02   0.167002E+01  -0.167002E+01  -0.177636E-14   0.381133E+00   0.110736E+01;
    0  0.286586E+02   0.168947E+01  -0.168947E+01  -0.384137E-13   0.379717E+00   0.113264E+01;
    0  0.289517E+02   0.170893E+01  -0.170893E+01   0.557332E-13   0.378287E+00   0.115809E+01;
    0  0.292411E+02   0.172862E+01  -0.172862E+01   0.424105E-13   0.376845E+00   0.118369E+01;
    4  0.290000E+02   0.171216E+01  -0.171216E+01  -0.295131E-10   0.378048E+00   0.116233E+01;
    

    % symm_re_40
    
    0  0.292888E+02   0.173197E+01  -0.173197E+01  -0.255465E-10   0.376604E+00   0.118796E+01;
    0  0.295739E+02   0.175214E+01  -0.175214E+01   0.224265E-12   0.375150E+00   0.121374E+01;
    0  0.298554E+02   0.177233E+01  -0.177233E+01   0.122347E-12   0.373685E+00   0.123965E+01;
    0  0.301331E+02   0.179253E+01  -0.179253E+01  -0.197620E-13   0.372212E+00   0.126568E+01;
    0  0.304071E+02   0.181276E+01  -0.181276E+01  -0.113243E-13   0.370731E+00   0.129182E+01;
    0  0.306774E+02   0.183300E+01  -0.183300E+01   0.222045E-14   0.369243E+00   0.131806E+01;
    0  0.309440E+02   0.185335E+01  -0.185335E+01   0.162093E-13   0.367750E+00   0.134439E+01;
    0  0.312069E+02   0.187442E+01  -0.187442E+01  -0.142109E-13   0.366253E+00   0.137081E+01;
    0  0.314660E+02   0.189565E+01  -0.189565E+01   0.821565E-14   0.364752E+00   0.139730E+01;
    0  0.317213E+02   0.191690E+01  -0.191690E+01   0.122125E-13   0.363248E+00   0.142387E+01;
    0  0.319730E+02   0.193817E+01  -0.193817E+01   0.146549E-13   0.361743E+00   0.145050E+01;
    0  0.322209E+02   0.195946E+01  -0.195946E+01  -0.355271E-14   0.360238E+00   0.147718E+01;
    0  0.324650E+02   0.198077E+01  -0.198077E+01   0.244249E-14   0.358733E+00   0.150392E+01;
    0  0.327055E+02   0.200216E+01  -0.200216E+01  -0.399680E-14   0.357229E+00   0.153071E+01;
    0  0.329422E+02   0.202415E+01  -0.202415E+01  -0.102141E-13   0.355728E+00   0.155755E+01;
    0  0.331751E+02   0.204614E+01  -0.204614E+01  -0.137668E-13   0.354230E+00   0.158443E+01;
    0  0.334044E+02   0.206816E+01  -0.206816E+01   0.177636E-14   0.352736E+00   0.161135E+01;
    0  0.336300E+02   0.209018E+01  -0.209018E+01  -0.799361E-14   0.351247E+00   0.163830E+01;
    0  0.338520E+02   0.211221E+01  -0.211221E+01  -0.102141E-13   0.349764E+00   0.166529E+01;
    0  0.340702E+02   0.213426E+01  -0.213426E+01  -0.266454E-14   0.348287E+00   0.169231E+01;
    0  0.342848E+02   0.215635E+01  -0.215635E+01  -0.754952E-14   0.346818E+00   0.171936E+01;
    0  0.344958E+02   0.217904E+01  -0.217904E+01  -0.222045E-14   0.345358E+00   0.174644E+01;
    0  0.347032E+02   0.220173E+01  -0.220173E+01   0.888178E-15   0.343906E+00   0.177354E+01;
    0  0.349070E+02   0.222442E+01  -0.222442E+01  -0.577316E-14   0.342464E+00   0.180067E+01;
    0  0.351072E+02   0.224712E+01  -0.224712E+01   0.000000E+00   0.341033E+00   0.182783E+01;
    0  0.353039E+02   0.226981E+01  -0.226981E+01   0.577316E-14   0.339613E+00   0.185500E+01;
    0  0.354971E+02   0.229250E+01  -0.229250E+01   0.355271E-14   0.338206E+00   0.188220E+01;
    0  0.356868E+02   0.231519E+01  -0.231519E+01  -0.310862E-14   0.336811E+00   0.190941E+01;
    0  0.358730E+02   0.233830E+01  -0.233830E+01  -0.399680E-14   0.335429E+00   0.193665E+01;
    0  0.360558E+02   0.236160E+01  -0.236160E+01   0.133227E-14   0.334062E+00   0.196390E+01;
    0  0.362353E+02   0.238489E+01  -0.238489E+01  -0.355271E-14   0.332710E+00   0.199116E+01;
    0  0.364113E+02   0.240817E+01  -0.240817E+01  -0.444089E-15   0.331373E+00   0.201845E+01;
    0  0.365840E+02   0.243144E+01  -0.243144E+01   0.266454E-14   0.330052E+00   0.204574E+01;
    0  0.367534E+02   0.245470E+01  -0.245470E+01   0.266454E-14   0.328748E+00   0.207305E+01;
    0  0.369196E+02   0.247794E+01  -0.247794E+01  -0.532907E-14   0.327462E+00   0.210037E+01;
    0  0.370825E+02   0.250116E+01  -0.250116E+01  -0.444089E-14   0.326193E+00   0.212770E+01;
    0  0.372422E+02   0.252495E+01  -0.252495E+01   0.888178E-15   0.324943E+00   0.215504E+01;
    0  0.373988E+02   0.254875E+01  -0.254875E+01   0.000000E+00   0.323711E+00   0.218239E+01;
    0  0.375522E+02   0.257254E+01  -0.257254E+01  -0.222045E-14   0.322500E+00   0.220974E+01;
    0  0.377026E+02   0.259630E+01  -0.259630E+01   0.106581E-13   0.321308E+00   0.223710E+01;
    0  0.378499E+02   0.262004E+01  -0.262004E+01  -0.177636E-14   0.320137E+00   0.226447E+01;
    0  0.379942E+02   0.264385E+01  -0.264385E+01  -0.666134E-14   0.318987E+00   0.229184E+01;
    0  0.381356E+02   0.266800E+01  -0.266800E+01   0.442757E-12   0.317858E+00   0.231920E+01;
    0  0.382740E+02   0.269214E+01  -0.269214E+01  -0.687894E-12   0.316752E+00   0.234657E+01;
    0  0.384095E+02   0.271651E+01  -0.271651E+01   0.112621E-11   0.315668E+00   0.237394E+01;
    0  0.385422E+02   0.274121E+01  -0.274121E+01  -0.912248E-11   0.314607E+00   0.240130E+01;
    0  0.386721E+02   0.276589E+01  -0.276589E+01  -0.180300E-12   0.313569E+00   0.242866E+01;
    0  0.387993E+02   0.279054E+01  -0.279054E+01  -0.140776E-12   0.312555E+00   0.245602E+01;
    0  0.389237E+02   0.281516E+01  -0.281516E+01   0.293987E-12   0.311565E+00   0.248337E+01;
    0  0.390454E+02   0.283976E+01  -0.283976E+01   0.109992E-10   0.310600E+00   0.251071E+01;
    0  0.391645E+02   0.286433E+01  -0.286433E+01  -0.578648E-12   0.309659E+00   0.253804E+01;
    0  0.392810E+02   0.288887E+01  -0.288887E+01   0.815348E-12   0.308744E+00   0.256536E+01;
    0  0.393949E+02   0.291338E+01  -0.291338E+01  -0.103788E-10   0.307854E+00   0.259267E+01;
    0  0.395064E+02   0.293840E+01  -0.293840E+01   0.233147E-12   0.306990E+00   0.261997E+01;
    0  0.396153E+02   0.296345E+01  -0.296345E+01  -0.166978E-10   0.306152E+00   0.264725E+01;
    0  0.397218E+02   0.298847E+01  -0.298847E+01   0.306422E-13   0.305341E+00   0.267452E+01;
    0  0.398260E+02   0.301346E+01  -0.301346E+01  -0.427214E-12   0.304557E+00   0.270177E+01;
    0  0.399278E+02   0.303841E+01  -0.303841E+01   0.399014E-11   0.303800E+00   0.272900E+01;
    0  0.400272E+02   0.306332E+01  -0.306332E+01  -0.608891E-11   0.303070E+00   0.275621E+01;
    4  0.400000E+02   0.305645E+01  -0.305645E+01  -0.347278E-12   0.303269E+00   0.274870E+01;
    ],

    % 07_re_40_par19_0.0
    [
    0  0.315426E+02   0.186278E+01  -0.156519E+01   0.297594E+00   0.413522E+00   0.959985E+00;
    0  0.316982E+02   0.186607E+01  -0.155727E+01   0.308793E+00   0.415499E+00   0.944339E+00;
    0  0.318602E+02   0.186973E+01  -0.154940E+01   0.320333E+00   0.417497E+00   0.928355E+00;
    0  0.320285E+02   0.187317E+01  -0.154244E+01   0.330722E+00   0.419514E+00   0.912051E+00;
    0  0.322032E+02   0.187639E+01  -0.153524E+01   0.341147E+00   0.421547E+00   0.895444E+00;
    0  0.323845E+02   0.187939E+01  -0.152779E+01   0.351599E+00   0.423594E+00   0.878554E+00;
    0  0.325723E+02   0.188220E+01  -0.152013E+01   0.362073E+00   0.425652E+00   0.861400E+00;
    0  0.327667E+02   0.188481E+01  -0.151225E+01   0.372560E+00   0.427721E+00   0.844000E+00;
    0  0.329678E+02   0.188744E+01  -0.150417E+01   0.383275E+00   0.429797E+00   0.826375E+00;
    0  0.331757E+02   0.189053E+01  -0.149729E+01   0.393239E+00   0.431879E+00   0.808545E+00;
    0  0.333903E+02   0.189430E+01  -0.149024E+01   0.404061E+00   0.433965E+00   0.790530E+00;
    0  0.336118E+02   0.189791E+01  -0.148301E+01   0.414907E+00   0.436053E+00   0.772350E+00;
    0  0.338402E+02   0.190138E+01  -0.147561E+01   0.425770E+00   0.438141E+00   0.754026E+00;
    0  0.340755E+02   0.190470E+01  -0.146806E+01   0.436640E+00   0.440227E+00   0.735579E+00;
    0  0.343178E+02   0.190789E+01  -0.146037E+01   0.447512E+00   0.442309E+00   0.717029E+00;
    0  0.345671E+02   0.191095E+01  -0.145335E+01   0.457602E+00   0.444387E+00   0.698397E+00;
    0  0.348235E+02   0.191389E+01  -0.144676E+01   0.467137E+00   0.446459E+00   0.679704E+00;
    0  0.350870E+02   0.191673E+01  -0.144005E+01   0.476678E+00   0.448523E+00   0.660970E+00;
    0  0.353576E+02   0.192018E+01  -0.143324E+01   0.486937E+00   0.450577E+00   0.642214E+00;
    0  0.356353E+02   0.192395E+01  -0.142635E+01   0.497605E+00   0.452622E+00   0.623458E+00;
    0  0.359201E+02   0.192764E+01  -0.141939E+01   0.508254E+00   0.454655E+00   0.604719E+00;
    0  0.362121E+02   0.193125E+01  -0.141237E+01   0.518874E+00   0.456676E+00   0.586019E+00;
    0  0.365113E+02   0.193478E+01  -0.140639E+01   0.528389E+00   0.458684E+00   0.567375E+00;
    0  0.368176E+02   0.193827E+01  -0.140054E+01   0.537733E+00   0.460678E+00   0.548805E+00;
    0  0.371311E+02   0.194217E+01  -0.139465E+01   0.547519E+00   0.462658E+00   0.530328E+00;
    0  0.374518E+02   0.194603E+01  -0.138875E+01   0.557280E+00   0.464622E+00   0.511960E+00;
    0  0.377795E+02   0.194984E+01  -0.138284E+01   0.567008E+00   0.466571E+00   0.493719E+00;
    0  0.381143E+02   0.195363E+01  -0.137694E+01   0.576696E+00   0.468505E+00   0.475619E+00;
    0  0.384562E+02   0.195752E+01  -0.137107E+01   0.586453E+00   0.470422E+00   0.457677E+00;
    0  0.388052E+02   0.196228E+01  -0.136598E+01   0.596299E+00   0.472323E+00   0.439907E+00;
    0  0.391611E+02   0.196703E+01  -0.136132E+01   0.605710E+00   0.474208E+00   0.422322E+00;
    0  0.395239E+02   0.197178E+01  -0.135670E+01   0.615074E+00   0.476078E+00   0.404937E+00;
    0  0.398936E+02   0.197652E+01  -0.135214E+01   0.624384E+00   0.477931E+00   0.387762E+00;
    0  0.402700E+02   0.198128E+01  -0.134765E+01   0.633633E+00   0.479769E+00   0.370810E+00;
    4  0.400000E+02   0.197788E+01  -0.135085E+01   0.627022E+00   0.478455E+00   0.382915E+00;
    0  0.410430E+02   0.199084E+01  -0.133892E+01   0.651920E+00   0.483399E+00   0.337618E+00;
    0  0.414394E+02   0.199566E+01  -0.133471E+01   0.660943E+00   0.485193E+00   0.321396E+00;
    0  0.418422E+02   0.200050E+01  -0.133129E+01   0.669211E+00   0.486973E+00   0.305435E+00;
    0  0.422513E+02   0.200539E+01  -0.132832E+01   0.677073E+00   0.488740E+00   0.289743E+00;
    0  0.426666E+02   0.201031E+01  -0.132545E+01   0.684860E+00   0.490495E+00   0.274327E+00;
    0  0.430880E+02   0.201590E+01  -0.132272E+01   0.693178E+00   0.492238E+00   0.259193E+00;
    0  0.435154E+02   0.202182E+01  -0.132012E+01   0.701699E+00   0.493972E+00   0.244346E+00;
    0  0.439487E+02   0.202780E+01  -0.131768E+01   0.710124E+00   0.495696E+00   0.229791E+00;
    0  0.443876E+02   0.203408E+01  -0.131539E+01   0.718692E+00   0.497412E+00   0.215533E+00;
    ],

    % 09_re_45_boven
    [
    0  0.411552E+02   0.199221E+01  -0.133771E+01   0.654496E+00   0.483911E+00   0.332975E+00;
    0  0.415533E+02   0.199703E+01  -0.133354E+01   0.663495E+00   0.485701E+00   0.316827E+00;
    0  0.419579E+02   0.200189E+01  -0.133044E+01   0.671455E+00   0.487477E+00   0.300942E+00;
    0  0.423688E+02   0.200679E+01  -0.132749E+01   0.679297E+00   0.489240E+00   0.285328E+00;
    0  0.427859E+02   0.201172E+01  -0.132466E+01   0.687061E+00   0.490992E+00   0.269991E+00;
    0  0.432091E+02   0.201758E+01  -0.132197E+01   0.695612E+00   0.492733E+00   0.254938E+00;
    0  0.436381E+02   0.202352E+01  -0.131941E+01   0.704106E+00   0.494463E+00   0.240174E+00;
    0  0.440730E+02   0.202951E+01  -0.131701E+01   0.712503E+00   0.496185E+00   0.225703E+00;
    0  0.445135E+02   0.203592E+01  -0.131477E+01   0.721155E+00   0.497898E+00   0.211530E+00;
    0  0.449595E+02   0.204244E+01  -0.131270E+01   0.729738E+00   0.499605E+00   0.197658E+00;
    ],

    % 10_re_50_boven
    [
    0  0.454108E+02   0.204902E+01  -0.131091E+01   0.738107E+00   0.501306E+00   0.184089E+00;
    0  0.458673E+02   0.205567E+01  -0.130959E+01   0.746080E+00   0.503003E+00   0.170827E+00;
    0  0.463288E+02   0.206240E+01  -0.130890E+01   0.753498E+00   0.504698E+00   0.157873E+00;
    0  0.467952E+02   0.206920E+01  -0.130840E+01   0.760807E+00   0.506391E+00   0.145231E+00;
    0  0.472662E+02   0.207609E+01  -0.130809E+01   0.768000E+00   0.508084E+00   0.132900E+00;
    0  0.477417E+02   0.208307E+01  -0.130799E+01   0.775071E+00   0.509779E+00   0.120884E+00;
    0  0.482214E+02   0.209013E+01  -0.130812E+01   0.782014E+00   0.511478E+00   0.109183E+00;
    0  0.487050E+02   0.209729E+01  -0.130847E+01   0.788822E+00   0.513183E+00   0.977973E-01;
    0  0.491924E+02   0.210455E+01  -0.130906E+01   0.795489E+00   0.514896E+00   0.867279E-01;
    0  0.496830E+02   0.211191E+01  -0.131004E+01   0.801867E+00   0.516618E+00   0.759728E-01;
    

    % 11_re_57_boven
    
    0  0.502755E+02   0.212180E+01  -0.131159E+01   0.810212E+00   0.518699E+00   0.634749E-01;
    0  0.508711E+02   0.213180E+01  -0.131345E+01   0.818349E+00   0.520800E+00   0.514060E-01;
    0  0.514686E+02   0.214186E+01  -0.131557E+01   0.826285E+00   0.522920E+00   0.397326E-01;
    0  0.520662E+02   0.215188E+01  -0.131784E+01   0.834044E+00   0.525059E+00   0.283965E-01;
    0  0.526622E+02   0.216176E+01  -0.132009E+01   0.841666E+00   0.527211E+00   0.173111E-01;
    0  0.532548E+02   0.217132E+01  -0.132259E+01   0.848731E+00   0.529366E+00   0.636764E-02;
    0  0.538426E+02   0.218043E+01  -0.132469E+01   0.855741E+00   0.531510E+00  -0.454208E-02;
    0  0.544250E+02   0.218895E+01  -0.132618E+01   0.862770E+00   0.533627E+00  -0.154962E-01;
    0  0.550022E+02   0.219683E+01  -0.132699E+01   0.869847E+00   0.535705E+00  -0.265249E-01;
    0  0.555750E+02   0.220411E+01  -0.132712E+01   0.876985E+00   0.537735E+00  -0.376134E-01;
    

    % 12_re_63_boven
    
    0  0.575671E+02   0.222579E+01  -0.132405E+01   0.901737E+00   0.544462E+00  -0.762404E-01;
    0  0.581377E+02   0.223178E+01  -0.132242E+01   0.909357E+00   0.546287E+00  -0.870247E-01;
    0  0.587107E+02   0.223785E+01  -0.132061E+01   0.917245E+00   0.548078E+00  -0.976323E-01;
    0  0.592870E+02   0.224378E+01  -0.131867E+01   0.925104E+00   0.549840E+00  -0.108039E+00;
    0  0.598670E+02   0.224958E+01  -0.131673E+01   0.932851E+00   0.551575E+00  -0.118228E+00;
    0  0.604511E+02   0.225574E+01  -0.131518E+01   0.940567E+00   0.553287E+00  -0.128186E+00;
    0  0.610396E+02   0.226188E+01  -0.131360E+01   0.948280E+00   0.554980E+00  -0.137903E+00;
    0  0.616326E+02   0.226798E+01  -0.131208E+01   0.955901E+00   0.556655E+00  -0.147375E+00;
    0  0.622305E+02   0.227408E+01  -0.131077E+01   0.963312E+00   0.558315E+00  -0.156598E+00;
    0  0.628331E+02   0.228019E+01  -0.130950E+01   0.970688E+00   0.559962E+00  -0.165570E+00;
    

    % 13_re_69_boven
    
    0  0.634407E+02   0.228632E+01  -0.130829E+01   0.978024E+00   0.561597E+00  -0.174291E+00;
    0  0.640532E+02   0.229247E+01  -0.130715E+01   0.985322E+00   0.563223E+00  -0.182764E+00;
    0  0.646706E+02   0.229866E+01  -0.130608E+01   0.992578E+00   0.564839E+00  -0.190989E+00;
    0  0.652930E+02   0.230489E+01  -0.130510E+01   0.999792E+00   0.566447E+00  -0.198971E+00;
    0  0.659202E+02   0.231117E+01  -0.130421E+01   0.100696E+01   0.568048E+00  -0.206712E+00;
    0  0.665524E+02   0.231751E+01  -0.130341E+01   0.101409E+01   0.569642E+00  -0.214218E+00;
    0  0.671894E+02   0.232389E+01  -0.130271E+01   0.102118E+01   0.571231E+00  -0.221491E+00;
    0  0.678311E+02   0.233034E+01  -0.130212E+01   0.102822E+01   0.572814E+00  -0.228537E+00;
    0  0.684776E+02   0.233684E+01  -0.130162E+01   0.103522E+01   0.574392E+00  -0.235360E+00;
    0  0.691288E+02   0.234341E+01  -0.130151E+01   0.104190E+01   0.575966E+00  -0.241966E+00;
    ],

    % 15_hopf_bifurcation
    [
    0  0.688841E+02   0.234094E+01  -0.130144E+01   0.103950E+01   0.575376E+00  -0.239514E+00;
    0  0.689656E+02   0.234176E+01  -0.130146E+01   0.104030E+01   0.575573E+00  -0.240335E+00;
    0  0.690472E+02   0.234259E+01  -0.130148E+01   0.104110E+01   0.575769E+00  -0.241152E+00;
    3  0.689781E+02   0.234189E+01  -0.130146E+01   0.104043E+01   0.575603E+00  -0.240460E+00;
    ],

    % 08_reynold_terug_re_30_beneden_par19_0.0
    [
    0  0.312499E+02   0.185695E+01  -0.158021E+01   0.276740E+00   0.409639E+00   0.990195E+00;
    0  0.311128E+02   0.185367E+01  -0.158730E+01   0.266368E+00   0.407737E+00   0.100473E+01;
    0  0.309817E+02   0.185012E+01  -0.159408E+01   0.256041E+00   0.405866E+00   0.101885E+01;
    0  0.308566E+02   0.184632E+01  -0.160179E+01   0.244534E+00   0.404026E+00   0.103256E+01;
    0  0.307375E+02   0.184226E+01  -0.160945E+01   0.232809E+00   0.402221E+00   0.104584E+01;
    0  0.306241E+02   0.183811E+01  -0.161682E+01   0.221293E+00   0.400452E+00   0.105867E+01;
    0  0.305166E+02   0.183518E+01  -0.162388E+01   0.211300E+00   0.398722E+00   0.107105E+01;
    0  0.304148E+02   0.183202E+01  -0.163063E+01   0.201392E+00   0.397033E+00   0.108295E+01;
    0  0.303187E+02   0.182857E+01  -0.163721E+01   0.191352E+00   0.395388E+00   0.109437E+01;
    0  0.302283E+02   0.182480E+01  -0.164483E+01   0.179971E+00   0.393787E+00   0.110529E+01;
    0  0.301434E+02   0.182074E+01  -0.165213E+01   0.168601E+00   0.392234E+00   0.111572E+01;
    0  0.300640E+02   0.181636E+01  -0.165911E+01   0.157245E+00   0.390730E+00   0.112563E+01;
    0  0.299901E+02   0.181183E+01  -0.166576E+01   0.146075E+00   0.389277E+00   0.113502E+01;
    0  0.299217E+02   0.180829E+01  -0.167206E+01   0.136223E+00   0.387877E+00   0.114388E+01;
    0  0.298586E+02   0.180441E+01  -0.167802E+01   0.126392E+00   0.386532E+00   0.115221E+01;
    0  0.298009E+02   0.180022E+01  -0.168515E+01   0.115069E+00   0.385244E+00   0.115999E+01;
    0  0.297486E+02   0.179569E+01  -0.169195E+01   0.103737E+00   0.384015E+00   0.116722E+01;
    0  0.297015E+02   0.179084E+01  -0.169842E+01   0.924198E-01   0.382845E+00   0.117390E+01;
    0  0.296596E+02   0.178565E+01  -0.170453E+01   0.811235E-01   0.381737E+00   0.118002E+01;
    0  0.296230E+02   0.178088E+01  -0.171029E+01   0.705947E-01   0.380692E+00   0.118558E+01;
    0  0.295916E+02   0.177652E+01  -0.171569E+01   0.608301E-01   0.379712E+00   0.119058E+01;
    0  0.295654E+02   0.177181E+01  -0.172202E+01   0.497925E-01   0.378797E+00   0.119501E+01;
    0  0.295444E+02   0.176678E+01  -0.172825E+01   0.385302E-01   0.377949E+00   0.119887E+01;
    0  0.295285E+02   0.176158E+01  -0.173412E+01   0.274573E-01   0.377169E+00   0.120216E+01;
    0  0.295177E+02   0.175603E+01  -0.173964E+01   0.163875E-01   0.376459E+00   0.120488E+01;
    0  0.295121E+02   0.175044E+01  -0.174480E+01   0.563376E-02   0.375818E+00   0.120703E+01;
    0  0.295117E+02   0.174571E+01  -0.174960E+01  -0.389224E-02   0.375247E+00   0.120861E+01;
    0  0.295163E+02   0.174061E+01  -0.175498E+01  -0.143645E-01   0.374749E+00   0.120963E+01;
    0  0.295261E+02   0.173516E+01  -0.176059E+01  -0.254331E-01   0.374322E+00   0.121008E+01;
    0  0.295411E+02   0.172935E+01  -0.176586E+01  -0.365058E-01   0.373968E+00   0.120997E+01;
    0  0.295612E+02   0.172318E+01  -0.177091E+01  -0.477301E-01   0.373686E+00   0.120930E+01;
    0  0.295865E+02   0.171667E+01  -0.177568E+01  -0.590138E-01   0.373478E+00   0.120809E+01;
    0  0.296169E+02   0.171130E+01  -0.178011E+01  -0.688080E-01   0.373343E+00   0.120632E+01;
    0  0.296525E+02   0.170561E+01  -0.178467E+01  -0.790608E-01   0.373281E+00   0.120401E+01;
    0  0.296934E+02   0.169956E+01  -0.178991E+01  -0.903529E-01   0.373292E+00   0.120117E+01;
    0  0.297396E+02   0.169316E+01  -0.179483E+01  -0.101666E+00   0.373376E+00   0.119780E+01;
    0  0.297910E+02   0.168642E+01  -0.179941E+01  -0.112997E+00   0.373532E+00   0.119390E+01;
    0  0.298477E+02   0.167933E+01  -0.180367E+01  -0.124338E+00   0.373761E+00   0.118950E+01;
    0  0.299098E+02   0.167318E+01  -0.180760E+01  -0.134424E+00   0.374060E+00   0.118459E+01;
    0  0.299772E+02   0.166694E+01  -0.181121E+01  -0.144273E+00   0.374431E+00   0.117918E+01;
    

    % 09_reynold_terug_re37.8_par19_0.0
    
    0  0.300501E+02   0.166035E+01  -0.181552E+01  -0.155171E+00   0.374871E+00   0.117329E+01;
    0  0.301285E+02   0.165343E+01  -0.181996E+01  -0.166523E+00   0.375380E+00   0.116692E+01;
    0  0.302123E+02   0.164619E+01  -0.182408E+01  -0.177892E+00   0.375956E+00   0.116010E+01;
    0  0.303018E+02   0.163863E+01  -0.182790E+01  -0.189271E+00   0.376600E+00   0.115282E+01;
    0  0.303968E+02   0.163183E+01  -0.183141E+01  -0.199583E+00   0.377308E+00   0.114510E+01;
    0  0.304976E+02   0.162514E+01  -0.183463E+01  -0.209487E+00   0.378079E+00   0.113695E+01;
    0  0.306041E+02   0.161813E+01  -0.183754E+01  -0.219412E+00   0.378913E+00   0.112839E+01;
    0  0.307163E+02   0.161082E+01  -0.184149E+01  -0.230668E+00   0.379807E+00   0.111943E+01;
    0  0.308344E+02   0.160321E+01  -0.184560E+01  -0.242389E+00   0.380759E+00   0.111009E+01;
    0  0.309584E+02   0.159532E+01  -0.184945E+01  -0.254132E+00   0.381767E+00   0.110037E+01;
    0  0.310884E+02   0.158856E+01  -0.185304E+01  -0.264476E+00   0.382830E+00   0.109030E+01;
    0  0.312244E+02   0.158153E+01  -0.185637E+01  -0.274840E+00   0.383945E+00   0.107989E+01;
    0  0.313665E+02   0.157421E+01  -0.185945E+01  -0.285244E+00   0.385109E+00   0.106915E+01;
    0  0.315148E+02   0.156661E+01  -0.186229E+01  -0.295682E+00   0.386321E+00   0.105811E+01;
    0  0.316693E+02   0.155874E+01  -0.186537E+01  -0.306633E+00   0.387578E+00   0.104677E+01;
    0  0.318301E+02   0.155064E+01  -0.186908E+01  -0.318438E+00   0.388878E+00   0.103516E+01;
    0  0.319972E+02   0.154373E+01  -0.187255E+01  -0.328820E+00   0.390217E+00   0.102330E+01;
    0  0.321708E+02   0.153657E+01  -0.187581E+01  -0.339239E+00   0.391592E+00   0.101119E+01;
    0  0.323508E+02   0.152917E+01  -0.187886E+01  -0.349686E+00   0.393002E+00   0.998862E+00;
    0  0.325374E+02   0.152154E+01  -0.188170E+01  -0.360157E+00   0.394444E+00   0.986330E+00;
    0  0.327307E+02   0.151370E+01  -0.188434E+01  -0.370642E+00   0.395913E+00   0.973612E+00;
    0  0.329305E+02   0.150566E+01  -0.188693E+01  -0.381272E+00   0.397408E+00   0.960727E+00;
    0  0.331372E+02   0.149856E+01  -0.188983E+01  -0.391264E+00   0.398925E+00   0.947693E+00;
    0  0.333506E+02   0.149154E+01  -0.189362E+01  -0.402081E+00   0.400461E+00   0.934527E+00;
    0  0.335708E+02   0.148434E+01  -0.189726E+01  -0.412923E+00   0.402013E+00   0.921248E+00;
    0  0.337979E+02   0.147697E+01  -0.190075E+01  -0.423783E+00   0.403579E+00   0.907875E+00;
    0  0.340320E+02   0.146945E+01  -0.190410E+01  -0.434653E+00   0.405154E+00   0.894425E+00;
    0  0.342730E+02   0.146179E+01  -0.190731E+01  -0.445525E+00   0.406735E+00   0.880917E+00;
    0  0.345210E+02   0.145454E+01  -0.191040E+01  -0.455859E+00   0.408321E+00   0.867367E+00;
    0  0.347761E+02   0.144797E+01  -0.191336E+01  -0.465393E+00   0.409906E+00   0.853794E+00;
    0  0.350383E+02   0.144128E+01  -0.191622E+01  -0.474933E+00   0.411489E+00   0.840215E+00;
    0  0.353076E+02   0.143449E+01  -0.191948E+01  -0.484985E+00   0.413066E+00   0.826646E+00;
    0  0.355840E+02   0.142761E+01  -0.192327E+01  -0.495656E+00   0.414635E+00   0.813104E+00;
    0  0.358675E+02   0.142066E+01  -0.192697E+01  -0.506309E+00   0.416191E+00   0.799605E+00;
    0  0.361582E+02   0.141366E+01  -0.193059E+01  -0.516935E+00   0.417733E+00   0.786165E+00;
    0  0.364561E+02   0.140746E+01  -0.193414E+01  -0.526683E+00   0.419257E+00   0.772799E+00;
    0  0.367611E+02   0.140161E+01  -0.193762E+01  -0.536010E+00   0.420761E+00   0.759522E+00;
    0  0.370733E+02   0.139573E+01  -0.194146E+01  -0.545731E+00   0.422241E+00   0.746348E+00;
    0  0.373926E+02   0.138983E+01  -0.194532E+01  -0.555498E+00   0.423695E+00   0.733289E+00;
    0  0.377190E+02   0.138392E+01  -0.194915E+01  -0.565232E+00   0.425120E+00   0.720361E+00;
    

    % 10_re40_par_19_0.0_onderaan
    
    0  0.380526E+02   0.137801E+01  -0.195294E+01  -0.574928E+00   0.426515E+00   0.707573E+00;
    0  0.383932E+02   0.137214E+01  -0.195671E+01  -0.584576E+00   0.427875E+00   0.694940E+00;
    0  0.387408E+02   0.136684E+01  -0.196141E+01  -0.594574E+00   0.429199E+00   0.682471E+00;
    0  0.390955E+02   0.136217E+01  -0.196616E+01  -0.603992E+00   0.430485E+00   0.670177E+00;
    0  0.394570E+02   0.135754E+01  -0.197091E+01  -0.613365E+00   0.431730E+00   0.658068E+00;
    0  0.398255E+02   0.135297E+01  -0.197565E+01  -0.622686E+00   0.432932E+00   0.646153E+00;
    0  0.402007E+02   0.134846E+01  -0.198041E+01  -0.631947E+00   0.434089E+00   0.634441E+00;
    4  0.400000E+02   0.135085E+01  -0.197788E+01  -0.627022E+00   0.433479E+00   0.640653E+00;
    ]
  };

  pitchfork = [
    3  0.295112E+02   0.174768E+01  -0.174768E+01   0.198612E-08   0.375472E+00   0.120803E+01
  ];

  hopf = [
    3  0.689781E+02   0.234189E+01  -0.130146E+01   0.104043E+01   0.575603E+00  -0.240460E+00
  ];


  legend_data = {};

  for i = 1 : length(data_sets)
    data = data_sets{i};
    plot(data(:, 2), data(:, column), 'b.-'); hold on;
    % interesting = find(data(:, 1) == 4 | data(:, 1) == 3);
    % scatter(data(interesting, 2), data(interesting, column), 'r'); hold on;
  end

  % Pitch
  pitch_f = scatter(pitchfork(2), pitchfork(column), 'r');

  % Hopf
  hopf_f = scatter(hopf(2), hopf(column), 'g');

  legend([pitch_f, hopf_f], {'Pitchfork','Hopf'});

  xlabel('Reynolds number')
  ylabel(names(column))
  grid on

  hold off;
end