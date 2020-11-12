import unittest
import Bio
from script import *

class TestSummarizeContents(unittest.TestCase):

    def test_summarize_contents(self):
        s = summarize_contents("/home/omar/ejercicio-biopython/data/AF323668.gbk")
        self.assertEqual(('name:AF323668.gbk', 'path:/home/omar/ejercicio-biopython/data/AF323668.gbk', 'num_records: 1 records', ['-id:AF323668.1', 'name:AF323668', 'Description:Bacteriophage bIL285, complete genome']), s)

        s = summarize_contents("/home/omar/ejercicio-biopython/data/ls_orchid.gbk")
        self.assertEqual(('name:ls_orchid.gbk', 'path:/home/omar/ejercicio-biopython/data/ls_orchid.gbk', 'num_records: 94 records', ['-id:Z78533.1', 'name:Z78533', 'Description:C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78532.1', 'name:Z78532', 'Description:C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78531.1', 'name:Z78531', 'Description:C.fasciculatum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78530.1', 'name:Z78530', 'Description:C.margaritaceum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78529.1', 'name:Z78529', 'Description:C.lichiangense 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78527.1', 'name:Z78527', 'Description:C.yatabeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78526.1', 'name:Z78526', 'Description:C.guttatum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78525.1', 'name:Z78525', 'Description:C.acaule 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78524.1', 'name:Z78524', 'Description:C.formosanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78523.1', 'name:Z78523', 'Description:C.himalaicum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78522.1', 'name:Z78522', 'Description:C.macranthum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78521.1', 'name:Z78521', 'Description:C.calceolus 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78520.1', 'name:Z78520', 'Description:C.segawai 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78519.1', 'name:Z78519', 'Description:C.pubescens 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78518.1', 'name:Z78518', 'Description:C.reginae 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78517.1', 'name:Z78517', 'Description:C.flavum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78516.1', 'name:Z78516', 'Description:C.passerinum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78515.1', 'name:Z78515', 'Description:M.xerophyticum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78514.1', 'name:Z78514', 'Description:P.schlimii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78513.1', 'name:Z78513', 'Description:P.besseae 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78512.1', 'name:Z78512', 'Description:P.wallisii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78511.1', 'name:Z78511', 'Description:P.exstaminodium 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78510.1', 'name:Z78510', 'Description:P.caricinum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78509.1', 'name:Z78509', 'Description:P.pearcei 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78508.1', 'name:Z78508', 'Description:P.longifolium 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78507.1', 'name:Z78507', 'Description:P.lindenii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78506.1', 'name:Z78506', 'Description:P.lindleyanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78505.1', 'name:Z78505', 'Description:P.sargentianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78504.1', 'name:Z78504', 'Description:P.kaiteurum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78503.1', 'name:Z78503', 'Description:P.czerwiakowianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78502.1', 'name:Z78502', 'Description:P.boissierianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78501.1', 'name:Z78501', 'Description:P.caudatum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78500.1', 'name:Z78500', 'Description:P.warszewiczianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78499.1', 'name:Z78499', 'Description:P.micranthum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78498.1', 'name:Z78498', 'Description:P.malipoense 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78497.1', 'name:Z78497', 'Description:P.delenatii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78496.1', 'name:Z78496', 'Description:P.armeniacum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78495.1', 'name:Z78495', 'Description:P.emersonii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78494.1', 'name:Z78494', 'Description:P.niveum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78493.1', 'name:Z78493', 'Description:P.godefroyae 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78492.1', 'name:Z78492', 'Description:P.bellatulum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78491.1', 'name:Z78491', 'Description:P.concolor 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78490.1', 'name:Z78490', 'Description:P.fairrieanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78489.1', 'name:Z78489', 'Description:P.druryi 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78488.1', 'name:Z78488', 'Description:P.tigrinum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78487.1', 'name:Z78487', 'Description:P.hirsutissimum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78486.1', 'name:Z78486', 'Description:P.barbigerum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78485.1', 'name:Z78485', 'Description:P.henryanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78484.1', 'name:Z78484', 'Description:P.charlesworthii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78483.1', 'name:Z78483', 'Description:P.villosum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78482.1', 'name:Z78482', 'Description:P.exul 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78481.1', 'name:Z78481', 'Description:P.insigne 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78480.1', 'name:Z78480', 'Description:P.gratrixianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78479.1', 'name:Z78479', 'Description:P.primulinum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78478.1', 'name:Z78478', 'Description:P.victoria 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78477.1', 'name:Z78477', 'Description:P.victoria 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78476.1', 'name:Z78476', 'Description:P.glaucophyllum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78475.1', 'name:Z78475', 'Description:P.supardii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78474.1', 'name:Z78474', 'Description:P.kolopakingii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78473.1', 'name:Z78473', 'Description:P.sanderianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78472.1', 'name:Z78472', 'Description:P.lowii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78471.1', 'name:Z78471', 'Description:P.dianthum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78470.1', 'name:Z78470', 'Description:P.parishii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78469.1', 'name:Z78469', 'Description:P.haynaldianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78468.1', 'name:Z78468', 'Description:P.adductum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78467.1', 'name:Z78467', 'Description:P.stonei 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78466.1', 'name:Z78466', 'Description:P.philippinense 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78465.1', 'name:Z78465', 'Description:P.rothschildianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78464.1', 'name:Z78464', 'Description:P.glanduliferum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78463.1', 'name:Z78463', 'Description:P.glanduliferum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78462.1', 'name:Z78462', 'Description:P.sukhakulii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78461.1', 'name:Z78461', 'Description:P.wardii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78460.1', 'name:Z78460', 'Description:P.ciliolare 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78459.1', 'name:Z78459', 'Description:P.dayanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78458.1', 'name:Z78458', 'Description:P.hennisianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78457.1', 'name:Z78457', 'Description:P.callosum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78456.1', 'name:Z78456', 'Description:P.tonsum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78455.1', 'name:Z78455', 'Description:P.javanicum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78454.1', 'name:Z78454', 'Description:P.fowliei 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78453.1', 'name:Z78453', 'Description:P.schoseri 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78452.1', 'name:Z78452', 'Description:P.bougainvilleanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78451.1', 'name:Z78451', 'Description:P.hookerae 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78450.1', 'name:Z78450', 'Description:P.papuanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78449.1', 'name:Z78449', 'Description:P.mastersianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78448.1', 'name:Z78448', 'Description:P.argus 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78447.1', 'name:Z78447', 'Description:P.venustum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78446.1', 'name:Z78446', 'Description:P.acmodontum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78445.1', 'name:Z78445', 'Description:P.urbanianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78444.1', 'name:Z78444', 'Description:P.appletonianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78443.1', 'name:Z78443', 'Description:P.lawrenceanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78442.1', 'name:Z78442', 'Description:P.bullenianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78441.1', 'name:Z78441', 'Description:P.superbiens 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78440.1', 'name:Z78440', 'Description:P.purpuratum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:Z78439.1', 'name:Z78439', 'Description:P.barbatum 5.8S rRNA gene and ITS1 and ITS2 DNA']),s)

        s = summarize_contents("/home/omar/ejercicio-biopython/data/NC_002703.gbk")
        self.assertEqual(('name:NC_002703.gbk', 'path:/home/omar/ejercicio-biopython/data/NC_002703.gbk', 'num_records: 1 records', ['-id:NC_002703.1', 'name:NC_002703', 'Description:Lactococcus phage Tuc2009, complete genome']),s)
        
        s = summarize_contents("/home/omar/ejercicio-biopython/data/ls_orchid.fasta")
        self.assertEqual(('name:ls_orchid.fasta', 'path:/home/omar/ejercicio-biopython/data/ls_orchid.fasta', 'num_records: 94 records', ['-id:gi|2765658|emb|Z78533.1|CIZ78533', 'name:gi|2765658|emb|Z78533.1|CIZ78533', 'Description:gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765657|emb|Z78532.1|CCZ78532', 'name:gi|2765657|emb|Z78532.1|CCZ78532', 'Description:gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765656|emb|Z78531.1|CFZ78531', 'name:gi|2765656|emb|Z78531.1|CFZ78531', 'Description:gi|2765656|emb|Z78531.1|CFZ78531 C.fasciculatum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765655|emb|Z78530.1|CMZ78530', 'name:gi|2765655|emb|Z78530.1|CMZ78530', 'Description:gi|2765655|emb|Z78530.1|CMZ78530 C.margaritaceum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765654|emb|Z78529.1|CLZ78529', 'name:gi|2765654|emb|Z78529.1|CLZ78529', 'Description:gi|2765654|emb|Z78529.1|CLZ78529 C.lichiangense 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765652|emb|Z78527.1|CYZ78527', 'name:gi|2765652|emb|Z78527.1|CYZ78527', 'Description:gi|2765652|emb|Z78527.1|CYZ78527 C.yatabeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765651|emb|Z78526.1|CGZ78526', 'name:gi|2765651|emb|Z78526.1|CGZ78526', 'Description:gi|2765651|emb|Z78526.1|CGZ78526 C.guttatum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765650|emb|Z78525.1|CAZ78525', 'name:gi|2765650|emb|Z78525.1|CAZ78525', 'Description:gi|2765650|emb|Z78525.1|CAZ78525 C.acaule 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765649|emb|Z78524.1|CFZ78524', 'name:gi|2765649|emb|Z78524.1|CFZ78524', 'Description:gi|2765649|emb|Z78524.1|CFZ78524 C.formosanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765648|emb|Z78523.1|CHZ78523', 'name:gi|2765648|emb|Z78523.1|CHZ78523', 'Description:gi|2765648|emb|Z78523.1|CHZ78523 C.himalaicum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765647|emb|Z78522.1|CMZ78522', 'name:gi|2765647|emb|Z78522.1|CMZ78522', 'Description:gi|2765647|emb|Z78522.1|CMZ78522 C.macranthum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765646|emb|Z78521.1|CCZ78521', 'name:gi|2765646|emb|Z78521.1|CCZ78521', 'Description:gi|2765646|emb|Z78521.1|CCZ78521 C.calceolus 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765645|emb|Z78520.1|CSZ78520', 'name:gi|2765645|emb|Z78520.1|CSZ78520', 'Description:gi|2765645|emb|Z78520.1|CSZ78520 C.segawai 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765644|emb|Z78519.1|CPZ78519', 'name:gi|2765644|emb|Z78519.1|CPZ78519', 'Description:gi|2765644|emb|Z78519.1|CPZ78519 C.pubescens 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765643|emb|Z78518.1|CRZ78518', 'name:gi|2765643|emb|Z78518.1|CRZ78518', 'Description:gi|2765643|emb|Z78518.1|CRZ78518 C.reginae 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765642|emb|Z78517.1|CFZ78517', 'name:gi|2765642|emb|Z78517.1|CFZ78517', 'Description:gi|2765642|emb|Z78517.1|CFZ78517 C.flavum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765641|emb|Z78516.1|CPZ78516', 'name:gi|2765641|emb|Z78516.1|CPZ78516', 'Description:gi|2765641|emb|Z78516.1|CPZ78516 C.passerinum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765640|emb|Z78515.1|MXZ78515', 'name:gi|2765640|emb|Z78515.1|MXZ78515', 'Description:gi|2765640|emb|Z78515.1|MXZ78515 M.xerophyticum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765639|emb|Z78514.1|PSZ78514', 'name:gi|2765639|emb|Z78514.1|PSZ78514', 'Description:gi|2765639|emb|Z78514.1|PSZ78514 P.schlimii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765638|emb|Z78513.1|PBZ78513', 'name:gi|2765638|emb|Z78513.1|PBZ78513', 'Description:gi|2765638|emb|Z78513.1|PBZ78513 P.besseae 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765637|emb|Z78512.1|PWZ78512', 'name:gi|2765637|emb|Z78512.1|PWZ78512', 'Description:gi|2765637|emb|Z78512.1|PWZ78512 P.wallisii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765636|emb|Z78511.1|PEZ78511', 'name:gi|2765636|emb|Z78511.1|PEZ78511', 'Description:gi|2765636|emb|Z78511.1|PEZ78511 P.exstaminodium 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765635|emb|Z78510.1|PCZ78510', 'name:gi|2765635|emb|Z78510.1|PCZ78510', 'Description:gi|2765635|emb|Z78510.1|PCZ78510 P.caricinum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765634|emb|Z78509.1|PPZ78509', 'name:gi|2765634|emb|Z78509.1|PPZ78509', 'Description:gi|2765634|emb|Z78509.1|PPZ78509 P.pearcei 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765633|emb|Z78508.1|PLZ78508', 'name:gi|2765633|emb|Z78508.1|PLZ78508', 'Description:gi|2765633|emb|Z78508.1|PLZ78508 P.longifolium 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765632|emb|Z78507.1|PLZ78507', 'name:gi|2765632|emb|Z78507.1|PLZ78507', 'Description:gi|2765632|emb|Z78507.1|PLZ78507 P.lindenii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765631|emb|Z78506.1|PLZ78506', 'name:gi|2765631|emb|Z78506.1|PLZ78506', 'Description:gi|2765631|emb|Z78506.1|PLZ78506 P.lindleyanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765630|emb|Z78505.1|PSZ78505', 'name:gi|2765630|emb|Z78505.1|PSZ78505', 'Description:gi|2765630|emb|Z78505.1|PSZ78505 P.sargentianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765629|emb|Z78504.1|PKZ78504', 'name:gi|2765629|emb|Z78504.1|PKZ78504', 'Description:gi|2765629|emb|Z78504.1|PKZ78504 P.kaiteurum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765628|emb|Z78503.1|PCZ78503', 'name:gi|2765628|emb|Z78503.1|PCZ78503', 'Description:gi|2765628|emb|Z78503.1|PCZ78503 P.czerwiakowianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765627|emb|Z78502.1|PBZ78502', 'name:gi|2765627|emb|Z78502.1|PBZ78502', 'Description:gi|2765627|emb|Z78502.1|PBZ78502 P.boissierianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765626|emb|Z78501.1|PCZ78501', 'name:gi|2765626|emb|Z78501.1|PCZ78501', 'Description:gi|2765626|emb|Z78501.1|PCZ78501 P.caudatum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765625|emb|Z78500.1|PWZ78500', 'name:gi|2765625|emb|Z78500.1|PWZ78500', 'Description:gi|2765625|emb|Z78500.1|PWZ78500 P.warszewiczianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765624|emb|Z78499.1|PMZ78499', 'name:gi|2765624|emb|Z78499.1|PMZ78499', 'Description:gi|2765624|emb|Z78499.1|PMZ78499 P.micranthum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765623|emb|Z78498.1|PMZ78498', 'name:gi|2765623|emb|Z78498.1|PMZ78498', 'Description:gi|2765623|emb|Z78498.1|PMZ78498 P.malipoense 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765622|emb|Z78497.1|PDZ78497', 'name:gi|2765622|emb|Z78497.1|PDZ78497', 'Description:gi|2765622|emb|Z78497.1|PDZ78497 P.delenatii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765621|emb|Z78496.1|PAZ78496', 'name:gi|2765621|emb|Z78496.1|PAZ78496', 'Description:gi|2765621|emb|Z78496.1|PAZ78496 P.armeniacum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765620|emb|Z78495.1|PEZ78495', 'name:gi|2765620|emb|Z78495.1|PEZ78495', 'Description:gi|2765620|emb|Z78495.1|PEZ78495 P.emersonii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765619|emb|Z78494.1|PNZ78494', 'name:gi|2765619|emb|Z78494.1|PNZ78494', 'Description:gi|2765619|emb|Z78494.1|PNZ78494 P.niveum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765618|emb|Z78493.1|PGZ78493', 'name:gi|2765618|emb|Z78493.1|PGZ78493', 'Description:gi|2765618|emb|Z78493.1|PGZ78493 P.godefroyae 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765617|emb|Z78492.1|PBZ78492', 'name:gi|2765617|emb|Z78492.1|PBZ78492', 'Description:gi|2765617|emb|Z78492.1|PBZ78492 P.bellatulum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765616|emb|Z78491.1|PCZ78491', 'name:gi|2765616|emb|Z78491.1|PCZ78491', 'Description:gi|2765616|emb|Z78491.1|PCZ78491 P.concolor 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765615|emb|Z78490.1|PFZ78490', 'name:gi|2765615|emb|Z78490.1|PFZ78490', 'Description:gi|2765615|emb|Z78490.1|PFZ78490 P.fairrieanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765614|emb|Z78489.1|PDZ78489', 'name:gi|2765614|emb|Z78489.1|PDZ78489', 'Description:gi|2765614|emb|Z78489.1|PDZ78489 P.druryi 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765613|emb|Z78488.1|PTZ78488', 'name:gi|2765613|emb|Z78488.1|PTZ78488', 'Description:gi|2765613|emb|Z78488.1|PTZ78488 P.tigrinum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765612|emb|Z78487.1|PHZ78487', 'name:gi|2765612|emb|Z78487.1|PHZ78487', 'Description:gi|2765612|emb|Z78487.1|PHZ78487 P.hirsutissimum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765611|emb|Z78486.1|PBZ78486', 'name:gi|2765611|emb|Z78486.1|PBZ78486', 'Description:gi|2765611|emb|Z78486.1|PBZ78486 P.barbigerum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765610|emb|Z78485.1|PHZ78485', 'name:gi|2765610|emb|Z78485.1|PHZ78485', 'Description:gi|2765610|emb|Z78485.1|PHZ78485 P.henryanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765609|emb|Z78484.1|PCZ78484', 'name:gi|2765609|emb|Z78484.1|PCZ78484', 'Description:gi|2765609|emb|Z78484.1|PCZ78484 P.charlesworthii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765608|emb|Z78483.1|PVZ78483', 'name:gi|2765608|emb|Z78483.1|PVZ78483', 'Description:gi|2765608|emb|Z78483.1|PVZ78483 P.villosum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765607|emb|Z78482.1|PEZ78482', 'name:gi|2765607|emb|Z78482.1|PEZ78482', 'Description:gi|2765607|emb|Z78482.1|PEZ78482 P.exul 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765606|emb|Z78481.1|PIZ78481', 'name:gi|2765606|emb|Z78481.1|PIZ78481', 'Description:gi|2765606|emb|Z78481.1|PIZ78481 P.insigne 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765605|emb|Z78480.1|PGZ78480', 'name:gi|2765605|emb|Z78480.1|PGZ78480', 'Description:gi|2765605|emb|Z78480.1|PGZ78480 P.gratrixianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765604|emb|Z78479.1|PPZ78479', 'name:gi|2765604|emb|Z78479.1|PPZ78479', 'Description:gi|2765604|emb|Z78479.1|PPZ78479 P.primulinum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765603|emb|Z78478.1|PVZ78478', 'name:gi|2765603|emb|Z78478.1|PVZ78478', 'Description:gi|2765603|emb|Z78478.1|PVZ78478 P.victoria 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765602|emb|Z78477.1|PVZ78477', 'name:gi|2765602|emb|Z78477.1|PVZ78477', 'Description:gi|2765602|emb|Z78477.1|PVZ78477 P.victoria 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765601|emb|Z78476.1|PGZ78476', 'name:gi|2765601|emb|Z78476.1|PGZ78476', 'Description:gi|2765601|emb|Z78476.1|PGZ78476 P.glaucophyllum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765600|emb|Z78475.1|PSZ78475', 'name:gi|2765600|emb|Z78475.1|PSZ78475', 'Description:gi|2765600|emb|Z78475.1|PSZ78475 P.supardii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765599|emb|Z78474.1|PKZ78474', 'name:gi|2765599|emb|Z78474.1|PKZ78474', 'Description:gi|2765599|emb|Z78474.1|PKZ78474 P.kolopakingii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765598|emb|Z78473.1|PSZ78473', 'name:gi|2765598|emb|Z78473.1|PSZ78473', 'Description:gi|2765598|emb|Z78473.1|PSZ78473 P.sanderianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765597|emb|Z78472.1|PLZ78472', 'name:gi|2765597|emb|Z78472.1|PLZ78472', 'Description:gi|2765597|emb|Z78472.1|PLZ78472 P.lowii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765596|emb|Z78471.1|PDZ78471', 'name:gi|2765596|emb|Z78471.1|PDZ78471', 'Description:gi|2765596|emb|Z78471.1|PDZ78471 P.dianthum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765595|emb|Z78470.1|PPZ78470', 'name:gi|2765595|emb|Z78470.1|PPZ78470', 'Description:gi|2765595|emb|Z78470.1|PPZ78470 P.parishii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765594|emb|Z78469.1|PHZ78469', 'name:gi|2765594|emb|Z78469.1|PHZ78469', 'Description:gi|2765594|emb|Z78469.1|PHZ78469 P.haynaldianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765593|emb|Z78468.1|PAZ78468', 'name:gi|2765593|emb|Z78468.1|PAZ78468', 'Description:gi|2765593|emb|Z78468.1|PAZ78468 P.adductum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765592|emb|Z78467.1|PSZ78467', 'name:gi|2765592|emb|Z78467.1|PSZ78467', 'Description:gi|2765592|emb|Z78467.1|PSZ78467 P.stonei 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765591|emb|Z78466.1|PPZ78466', 'name:gi|2765591|emb|Z78466.1|PPZ78466', 'Description:gi|2765591|emb|Z78466.1|PPZ78466 P.philippinense 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765590|emb|Z78465.1|PRZ78465', 'name:gi|2765590|emb|Z78465.1|PRZ78465', 'Description:gi|2765590|emb|Z78465.1|PRZ78465 P.rothschildianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765589|emb|Z78464.1|PGZ78464', 'name:gi|2765589|emb|Z78464.1|PGZ78464', 'Description:gi|2765589|emb|Z78464.1|PGZ78464 P.glanduliferum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765588|emb|Z78463.1|PGZ78463', 'name:gi|2765588|emb|Z78463.1|PGZ78463', 'Description:gi|2765588|emb|Z78463.1|PGZ78463 P.glanduliferum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765587|emb|Z78462.1|PSZ78462', 'name:gi|2765587|emb|Z78462.1|PSZ78462', 'Description:gi|2765587|emb|Z78462.1|PSZ78462 P.sukhakulii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765586|emb|Z78461.1|PWZ78461', 'name:gi|2765586|emb|Z78461.1|PWZ78461', 'Description:gi|2765586|emb|Z78461.1|PWZ78461 P.wardii 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765585|emb|Z78460.1|PCZ78460', 'name:gi|2765585|emb|Z78460.1|PCZ78460', 'Description:gi|2765585|emb|Z78460.1|PCZ78460 P.ciliolare 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765584|emb|Z78459.1|PDZ78459', 'name:gi|2765584|emb|Z78459.1|PDZ78459', 'Description:gi|2765584|emb|Z78459.1|PDZ78459 P.dayanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765583|emb|Z78458.1|PHZ78458', 'name:gi|2765583|emb|Z78458.1|PHZ78458', 'Description:gi|2765583|emb|Z78458.1|PHZ78458 P.hennisianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765582|emb|Z78457.1|PCZ78457', 'name:gi|2765582|emb|Z78457.1|PCZ78457', 'Description:gi|2765582|emb|Z78457.1|PCZ78457 P.callosum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765581|emb|Z78456.1|PTZ78456', 'name:gi|2765581|emb|Z78456.1|PTZ78456', 'Description:gi|2765581|emb|Z78456.1|PTZ78456 P.tonsum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765580|emb|Z78455.1|PJZ78455', 'name:gi|2765580|emb|Z78455.1|PJZ78455', 'Description:gi|2765580|emb|Z78455.1|PJZ78455 P.javanicum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765579|emb|Z78454.1|PFZ78454', 'name:gi|2765579|emb|Z78454.1|PFZ78454', 'Description:gi|2765579|emb|Z78454.1|PFZ78454 P.fowliei 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765578|emb|Z78453.1|PSZ78453', 'name:gi|2765578|emb|Z78453.1|PSZ78453', 'Description:gi|2765578|emb|Z78453.1|PSZ78453 P.schoseri 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765577|emb|Z78452.1|PBZ78452', 'name:gi|2765577|emb|Z78452.1|PBZ78452', 'Description:gi|2765577|emb|Z78452.1|PBZ78452 P.bougainvilleanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765576|emb|Z78451.1|PHZ78451', 'name:gi|2765576|emb|Z78451.1|PHZ78451', 'Description:gi|2765576|emb|Z78451.1|PHZ78451 P.hookerae 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765575|emb|Z78450.1|PPZ78450', 'name:gi|2765575|emb|Z78450.1|PPZ78450', 'Description:gi|2765575|emb|Z78450.1|PPZ78450 P.papuanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765574|emb|Z78449.1|PMZ78449', 'name:gi|2765574|emb|Z78449.1|PMZ78449', 'Description:gi|2765574|emb|Z78449.1|PMZ78449 P.mastersianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765573|emb|Z78448.1|PAZ78448', 'name:gi|2765573|emb|Z78448.1|PAZ78448', 'Description:gi|2765573|emb|Z78448.1|PAZ78448 P.argus 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765572|emb|Z78447.1|PVZ78447', 'name:gi|2765572|emb|Z78447.1|PVZ78447', 'Description:gi|2765572|emb|Z78447.1|PVZ78447 P.venustum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765571|emb|Z78446.1|PAZ78446', 'name:gi|2765571|emb|Z78446.1|PAZ78446', 'Description:gi|2765571|emb|Z78446.1|PAZ78446 P.acmodontum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765570|emb|Z78445.1|PUZ78445', 'name:gi|2765570|emb|Z78445.1|PUZ78445', 'Description:gi|2765570|emb|Z78445.1|PUZ78445 P.urbanianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765569|emb|Z78444.1|PAZ78444', 'name:gi|2765569|emb|Z78444.1|PAZ78444', 'Description:gi|2765569|emb|Z78444.1|PAZ78444 P.appletonianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765568|emb|Z78443.1|PLZ78443', 'name:gi|2765568|emb|Z78443.1|PLZ78443', 'Description:gi|2765568|emb|Z78443.1|PLZ78443 P.lawrenceanum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765567|emb|Z78442.1|PBZ78442', 'name:gi|2765567|emb|Z78442.1|PBZ78442', 'Description:gi|2765567|emb|Z78442.1|PBZ78442 P.bullenianum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765566|emb|Z78441.1|PSZ78441', 'name:gi|2765566|emb|Z78441.1|PSZ78441', 'Description:gi|2765566|emb|Z78441.1|PSZ78441 P.superbiens 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765565|emb|Z78440.1|PPZ78440', 'name:gi|2765565|emb|Z78440.1|PPZ78440', 'Description:gi|2765565|emb|Z78440.1|PPZ78440 P.purpuratum 5.8S rRNA gene and ITS1 and ITS2 DNA', '-id:gi|2765564|emb|Z78439.1|PBZ78439', 'name:gi|2765564|emb|Z78439.1|PBZ78439', 'Description:gi|2765564|emb|Z78439.1|PBZ78439 P.barbatum 5.8S rRNA gene and ITS1 and ITS2 DNA']),s)

        s = summarize_contents("/home/omar/ejercicio-biopython/data/m_cold.fasta")
        self.assertEqual(('name:m_cold.fasta', 'path:/home/omar/ejercicio-biopython/data/m_cold.fasta', 'num_records: 1 records', ['-id:gi|8332116|gb|BE037100.1|BE037100', 'name:gi|8332116|gb|BE037100.1|BE037100', "Description:gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence"]),s)

        s = summarize_contents("/home/omar/ejercicio-biopython/data/opuntia.fasta")
        self.assertEqual(('name:opuntia.fasta', 'path:/home/omar/ejercicio-biopython/data/opuntia.fasta', 'num_records: 7 records', ['-id:gi|6273291|gb|AF191665.1|AF191665', 'name:gi|6273291|gb|AF191665.1|AF191665', 'Description:gi|6273291|gb|AF191665.1|AF191665 Opuntia marenae rpl16 gene; chloroplast gene for chloroplast product, partial intron sequence', '-id:gi|6273290|gb|AF191664.1|AF191664', 'name:gi|6273290|gb|AF191664.1|AF191664', 'Description:gi|6273290|gb|AF191664.1|AF191664 Opuntia clavata rpl16 gene; chloroplast gene for chloroplast product, partial intron sequence', '-id:gi|6273289|gb|AF191663.1|AF191663', 'name:gi|6273289|gb|AF191663.1|AF191663', 'Description:gi|6273289|gb|AF191663.1|AF191663 Opuntia bradtiana rpl16 gene; chloroplast gene for chloroplast product, partial intron sequence', '-id:gi|6273287|gb|AF191661.1|AF191661', 'name:gi|6273287|gb|AF191661.1|AF191661', 'Description:gi|6273287|gb|AF191661.1|AF191661 Opuntia kuehnrichiana rpl16 gene; chloroplast gene for chloroplast product, partial intron sequence', '-id:gi|6273286|gb|AF191660.1|AF191660', 'name:gi|6273286|gb|AF191660.1|AF191660', 'Description:gi|6273286|gb|AF191660.1|AF191660 Opuntia echinacea rpl16 gene; chloroplast gene for chloroplast product, partial intron sequence', '-id:gi|6273285|gb|AF191659.1|AF191659', 'name:gi|6273285|gb|AF191659.1|AF191659', 'Description:gi|6273285|gb|AF191659.1|AF191659 Opuntia pachypus rpl16 gene; chloroplast gene for chloroplast product, partial intron sequence', '-id:gi|6273284|gb|AF191658.1|AF191658', 'name:gi|6273284|gb|AF191658.1|AF191658', 'Description:gi|6273284|gb|AF191658.1|AF191658 Opuntia subulata rpl16 gene; chloroplast gene for chloroplast product, partial intron sequence']),s)

    #Función test para probar la función concatenate_and_get_reverse_of_complement de script.py
    def test_concatenate_and_get_reverse_of_complement(self):
        
        #Al parecer, la función del complemento inverso de biopython funciona bien para cualquier contenido en las cadenas
        actual = concatenate_and_get_reverse_of_complement(DNA_sequence_1, DNA_sequence_2)
        self.assertEqual("GCGATTTCGATCCTATATAGGCCCATCGATC", actual)

        actual = concatenate_and_get_reverse_of_complement("Hola", "mundo")
        self.assertEqual("ohnakuloD",actual)

        actual = concatenate_and_get_reverse_of_complement("@#$%", "^&*()")
        self.assertEqual(")(*&^%$#@", actual)

        actual = concatenate_and_get_reverse_of_complement("TTgCaAaACtgcACC", "GGCatTtaCcAAT")
        self.assertEqual("ATTgGtaAatGCCGGTgcaGTtTtGcAA", actual)
        
        actual = concatenate_and_get_reverse_of_complement("aa t cccg gactt", "gct ac tgg a")
        self.assertEqual("t cca gt agcaagtc cggg a tt", actual)

        actual = concatenate_and_get_reverse_of_complement("", "")
        self.assertEqual("", actual)

        actual = concatenate_and_get_reverse_of_complement("1234","5678")
        self.assertEqual("87654321", actual)

        self.assertRaises(TypeError, concatenate_and_get_reverse_of_complement,1234, 5678)

        self.assertRaises(TypeError, concatenate_and_get_reverse_of_complement, None, None)

    #Función para probar la función que imprime mRNA, proteínas y codones de paro, usando código estándar
    def test_print_protein_and_stop_codon_using_standard_table(self):

        actual = print_protein_and_stop_codon_using_standard_table("")
        self.assertEqual({'mRNA': Seq('')}, actual)

        actual = print_protein_and_stop_codon_using_standard_table("atttagcaggtttacgaccca")
        self.assertEqual({'mRNA': Seq('AUUUAGCAGGUUUACGACCCA')}, actual)

        actual = print_protein_and_stop_codon_using_standard_table("atgatttagcaggtttacgaccca")
        self.assertEqual({'mRNA': Seq('AUGAUUUAGCAGGUUUACGACCCA'), 'proteins': Seq('MI*QVYDP'), 'stop_codons': ['TAG']}, actual)

        self.assertRaises(Bio.Data.CodonTable.TranslationError, print_protein_and_stop_codon_using_standard_table, "GTGAAA AAGA TGCAATC TATCGT ACTCGCA CTT T C CCT")

        self.assertRaises(Bio.Data.CodonTable.TranslationError, print_protein_and_stop_codon_using_standard_table, "@#$%^^*()(*&^%$")

        self.assertRaises(TypeError, print_protein_and_stop_codon_using_standard_table, None)

        self.assertRaises(TypeError, print_protein_and_stop_codon_using_standard_table, 123456789)

        self.assertRaises(Bio.Data.CodonTable.TranslationError, print_protein_and_stop_codon_using_standard_table, "123456789")


        #Función para probar la función que imprime mRNA, proteínas y codones de paro, usando código estándar
        def test_print_proteins_and_codons_using_mitocondrial_yeast_table(self):

            actual = print_protein_and_stop_codon_using_standard_table("")
            self.assertEqual({'mRNA': Seq('')}, actual)

            actual = print_protein_and_stop_codon_using_standard_table("atttagcaggtttacgaccca")
            self.assertEqual({'mRNA': Seq('AUUUAGCAGGUUUACGACCCA')}, actual)

            actual = print_protein_and_stop_codon_using_standard_table("atgatttagcaggtttacgaccca")
            self.assertEqual({'mRNA': Seq('AUGAUUUAGCAGGUUUACGACCCA'), 'proteins': Seq('MI*QVYDP'), 'stop_codons': ['TAG']}, actual)

            self.assertRaises(Bio.Data.CodonTable.TranslationError, print_protein_and_stop_codon_using_standard_table, "GTGAAA AAGA TGCAATC TATCGT ACTCGCA CTT T C CCT")

            self.assertRaises(Bio.Data.CodonTable.TranslationError, print_protein_and_stop_codon_using_standard_table, "@#$%^^*()(*&^%$")

            self.assertRaises(TypeError, print_protein_and_stop_codon_using_standard_table, None)

            self.assertRaises(TypeError, print_protein_and_stop_codon_using_standard_table, 123456789)

            self.assertRaises(Bio.Data.CodonTable.TranslationError, print_protein_and_stop_codon_using_standard_table, "123456789")


if __name__ == "__main__":
    unittest.main()