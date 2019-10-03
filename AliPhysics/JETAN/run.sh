aliroot -b << EOF 
.x testua1.C("et14")
.q
EOF
mv jets.root jets_pt_0_r_1.0/150-180GeV.root


aliroot -b << EOF 
.x testua1.C("et13")
.q
EOF
mv jets.root jets_pt_0_r_1.0/125-150GeV.root

aliroot -b << EOF 
.x testua1.C("et12")
.q
EOF
mv jets.root jets_pt_0_r_1.0/104-125GeV.root

aliroot -b << EOF 
.x testua1.C("et11")
.q
EOF
mv jets.root jets_pt_0_r_1.0/86-104GeV.root


aliroot -b << EOF 
.x testua1.C("et10")
.q
EOF
mv jets.root  jets_pt_0_r_1.0/72-86GeV.root

aliroot -b << EOF 
.x testua1.C("et9")
.q
EOF
mv jets.root jets_pt_0_r_1.0/60-72GeV.root

aliroot -b << EOF 
.x testua1.C("et8")
.q
EOF
mv jets.root jets_pt_0_r_1.0/50-60GeV.root


aliroot -b << EOF 
.x testua1.C("et7")
.q
EOF
mv jets.root jets_pt_0_r_1.0/42-50GeV.root

aliroot -b << EOF 
.x testua1.C("et6")
.q
EOF
mv jets.root jets_pt_0_r_1.0/35-42GeV.root


aliroot -b << EOF 
.x testua1.C("et5")
.q
EOF
mv jets.root jets_pt_0_r_1.0/29-35GeV.root


aliroot -b << EOF 
.x testua1.C("et4")
.q
EOF
mv jets.root jets_pt_0_r_1.0/24-29GeV.root


aliroot -b << EOF 
.x testua1.C("et3")
.q
EOF
mv jets.root jets_pt_0_r_1.0/20-24GeV.root


