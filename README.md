# personalization2020
Dynamic binaural synthesis with mixed structural modeling: the personalization. /n
The personalization is a combination of HRTF selection (from a database, in this case, CIPIC) and a spherical head model with ear displacement. The different theories used in this personalization are described in https://www.dei.unipd.it/~geronazzo/data/publ/2020_ICASSP_camera_ready.pdf and https://projekter.aau.dk/projekter/files/291096949/Jason_Tissieres_SMC_2018.pdf. /n
Procedure:
1)	Download CIPIC database: https://www.ece.ucdavis.edu/cipic/spatial-sound/hrtf-data/.
2)	Select HRTF dataset with https://github.com/msmhrtf/sel.
3)	Run transformToTextWithDelaunay.m after entering HRTF dataset ID (CIPIC ID) and anthopometric dimensions (Head width X1 and Head depth X3)
4)	Output: personalised dataset in SOFA format and txt format. The txt file is formatted to be used in https://github.com/msmhrtf/binsdn.
