#!/bin/bash
for ((i=0; i<5; i++)); do ./build/MixedIntegerGOMEA_O -n100 -P1004 -Ooutput/meeting_25/rn_ew_1n/m3i1 -Irandomnetwork_ew_1nodes -F2_0_0_-1_0 -N1 -b0.99 -r0 -M150 -Y; done;
for ((i=0; i<5; i++)); do ./build/MixedIntegerGOMEA_O -n100 -P1004 -Ooutput/meeting_25/rn_ew_1n/m3i3 -Irandomnetwork_ew_1nodes -F2_0_0_-1_0 -N1 -b0.99 -r0 -M150 -Y -g; done;
for ((i=0; i<5; i++)); do ./build/MixedIntegerGOMEA_O -n100 -P1004 -Ooutput/meeting_25/rn_ew_1n/m3stari1 -Irandomnetwork_ew_1nodes -F2_0_0_-1_0 -N1 -b0.99 -r0 -M150 -Y -f; done;
for ((i=0; i<5; i++)); do ./build/MixedIntegerGOMEA_O -n100 -P1004 -Ooutput/meeting_25/rn_ew_1n/m3stari3 -Irandomnetwork_ew_1nodes -F2_0_0_-1_0 -N1 -b0.99 -r0 -M150 -Y -f -g; done;
for ((i=0; i<5; i++)); do ./build/MixedIntegerGOMEA_O -n100 -P1004 -Ooutput/meeting_25/n7/m3i1 -Inetwork7 -F2_0_0_-1_0 -N1 -b0.99 -r0 -M30 -Y; done;
for ((i=0; i<5; i++)); do ./build/MixedIntegerGOMEA_O -n100 -P1004 -Ooutput/meeting_25/n7/m3i3 -Inetwork7 -F2_0_0_-1_0 -N1 -b0.99 -r0 -M30 -Y -g; done;
for ((i=0; i<5; i++)); do ./build/MixedIntegerGOMEA_O -n100 -P1004 -Ooutput/meeting_25/n7/m3stari1 -Inetwork7 -F2_0_0_-1_0 -N1 -b0.99 -r0 -M30 -Y -f; done;
for ((i=0; i<5; i++)); do ./build/MixedIntegerGOMEA_O -n100 -P1004 -Ooutput/meeting_25/n7/m3stari3 -Inetwork7 -F2_0_0_-1_0 -N1 -b0.99 -r0 -M30 -Y -f -g; done;