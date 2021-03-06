#!/bin/bash
# : $(./pathtracer -t 8 -s 256 -l 4 -m 7 -c CBdragon_microfacet_au_05_cam_settings_head.txt -b 0.125000 -d 2.6 -r 480 360 -f images/4_dragon_f_01.png ../dae/sky/CBdragon_microfacet_au_05.dae >/dev/null 2>&1)
# : $(./pathtracer -t 8 -s 256 -l 4 -m 7 -c CBdragon_microfacet_au_05_cam_settings_body.txt -b 0.125000 -d 3.0 -r 480 360 -f images/4_dragon_f_02.png ../dae/sky/CBdragon_microfacet_au_05.dae >/dev/null 2>&1)
# : $(./pathtracer -t 8 -s 256 -l 4 -m 7 -c CBdragon_microfacet_au_05_cam_tail.txt -b 0.125000 -d 3.3 -r 480 360 -f images/4_dragon_f_03.png ../dae/sky/CBdragon_microfacet_au_05.dae >/dev/null 2>&1)
# : $(./pathtracer -t 8 -s 256 -l 4 -m 7 -c CBdragon_microfacet_au_05_cam_backwall.txt -b 0.125000 -d 4.2 -r 480 360 -f images/4_dragon_f_04.png ../dae/sky/CBdragon_microfacet_au_05.dae >/dev/null 2>&1)
# : $(./pathtracer -t 8 -s 256 -l 4 -m 7 -c CBdragon_microfacet_au_05_cam_settings_head.txt -b 0 -d 2.6 -r 480 360 -f images/4_dragon_a_01.png ../dae/sky/CBdragon_microfacet_au_05.dae >/dev/null 2>&1)
# : $(./pathtracer -t 8 -s 256 -l 4 -m 7 -c CBdragon_microfacet_au_05_cam_settings_head.txt -b 0.062500 -d 2.6 -r 480 360 -f images/4_dragon_a_02.png ../dae/sky/CBdragon_microfacet_au_05.dae >/dev/null 2>&1)
# : $(./pathtracer -t 8 -s 256 -l 4 -m 7 -c CBdragon_microfacet_au_05_cam_settings_head.txt -b 0.125000 -d 2.6 -r 480 360 -f images/4_dragon_a_03.png ../dae/sky/CBdragon_microfacet_au_05.dae >/dev/null 2>&1)
: $(./pathtracer -t 8 -s 1024 -l 4 -m 10 -r 480 360 -f cloud_a10_s09.png ../dae/sky/CBcloud.dae >/dev/null 2>&1)
: $(./pathtracer -t 8 -s 1024 -l 4 -m 10 -r 480 360 -f cloud_a10_s90.png ../dae/sky/CBcloud_s90.dae >/dev/null 2>&1)
