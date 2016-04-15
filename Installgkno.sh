if [ -e ./gkno_launcher/gkno ]
then
        echo "know already built: skipping"
else
        git clone https://github.com/gkno/gkno_launcher.git
        cd gkno_launcher/
        ./gkno build
        cd ..
fi

if [ -e ./gkno_launcher/resources/homo_sapiens/build_37_version_3/human_reference_v37_decoys.fa ]
then
        echo "human reference paramaters set already downloaded: skipping"
else
        cd gkno_launcher
        ./gkno add-resource human
        ./gkno bwa-index -r ./resources/homo_sapiens/build_37_version_3/human_reference_v37_decoys.fa -x ./resources/homo_sapiens/build_37_version_3/human_reference_v37_decoys
        cd ..
fi

cd ../../../
