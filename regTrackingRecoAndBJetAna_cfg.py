import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('standard')
ivars.register('jobIndex',mult=ivars.multiplicity.singleton,info="for testing")

ivars.jobIndex=99
ivars.files='file:/mnt/hadoop/cms/store/user/yenjie/MC_Production/Pythia80_HydjetDrum_mix01/RECO/set2_random10000_HydjetDrum_452.root'
ivars.output = 'test.root'

ivars.maxEvents = -1

ivars.parseArguments()

import FWCore.ParameterSet.Config as cms

isMC = True
hiReco = True
reReco = True
hasSimInfo = False
genTag = "hiSignal"
#hltFilter = "HLT_Jet40_v1"
#hltFilter = "HLT_HIL2Mu*_v*"
hltFilter = "HLT_HIJet80_v*"
trigResults = 'TriggerResults::RECO'
#trigResults = 'TriggerResults::HISIGNAL'
gTag = 'STARTHI44_V7::All'
#gTag = 'GR_R_44_V10::All'
hiMode = True
redoPFJets = False

# some important triggers:  HLT_Jet40_v1, HLT_HIL2Mu7_v1'

if hiReco:
    svTracks = "hiSecondaryVertexSelectedTracks"
    #svTracks = "hiGeneralTracks"
    pvProducer = "offlinePrimaryVertices"
    #print "hacked to look at hiGeneralTracks"
    #svTracks = "hiGeneralTracks"
    #pvProducer = "hiSelectedVertex"
else:
    svTracks = "generalTracks"
    pvProducer = "offlinePrimaryVertices"
    
print "Reco'ing SV's w/ ", svTracks, ", PV w/ ", pvProducer 

process = cms.Process('BJET')

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
#process.Timing = cms.Service("Timing")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')

if hiReco:
    process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
    process.load('Configuration.EventContent.EventContentHeavyIons_cff')
else:
    process.load('Configuration.StandardSequences.Reconstruction_cff')
    process.load('Configuration.EventContent.EventContent_cff')
    
if isMC:
    process.load('SimGeneral.MixingModule.mixNoPU_cfi')
    process.load('Configuration.StandardSequences.RawToDigi_cff')
else:
    process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
    
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(
    ivars.maxEvents
    #-1
    )
)

# Input source
process.source = cms.Source("PoolSource",
                            #secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring(
    ivars.files
    ),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            #eventsToProcess = cms.untracked.VEventRange('182838:7396593-182838:7396593'),
                            #skipEvents = cms.untracked.uint32(33)
                            )

process.options = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True)
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = gTag
print "global Tag = ", process.GlobalTag.globaltag

print "customizing JP calib "


process.jpCalib = cms.ESSource("PoolDBESSource",
                               DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
                               timetype = cms.string('runnumber'),
                               toGet = cms.VPSet(
    cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
             tag = cms.string('probBTagPDF2D_tag_mc'),
             #label = cms.untracked.string("IC5Calo")
             ),
    
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag    = cms.string('probBTagPDF3D_tag_mc'),
             #label = cms.untracked.string("AK3PF")
             ),
    
    ),
                               connect = cms.string("sqlite_file:jpCalib_HIMC44X_2011_v1.db"),
                               )
process.es_prefer_jpCalib = cms.ESPrefer('PoolDBESSource','jpCalib')


######################
# Hi specific reco
######################
# load centrality
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

if isMC:
    process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFtowers"),
	nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
	centralitySrc = cms.InputTag("hiCentrality")
	)
else:
    process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFtowers"),
	nonDefaultGlauberModel = cms.string(""),
	centralitySrc = cms.InputTag("hiCentrality")
	)

# need to produce centrality for other analyzers
process.load("RecoHI.HiCentralityAlgos.HiCentrality_cfi")
process.hiCentrality.producePixelTracks = cms.bool(False)
process.hiCentrality.produceETmidRapidity = cms.bool(False)
process.hiCentrality.producePixelhits = cms.bool(False)
process.hiCentrality.produceEcalhits = cms.bool(False)
process.hiCentrality.produceZDChits = cms.bool(False)
process.hiCentrality.produceBasicClusters = cms.bool(False)
process.hiCentrality.produceHFhits = cms.bool(True)
process.hiCentrality.produceTracks = cms.bool(False)


#if isMC==False:
#    
#    import PhysicsTools.PythonAnalysis.LumiList as LumiList
#    import FWCore.ParameterSet.Types as CfgTypes
#    myLumis = LumiList.LumiList(filename = 'json.txt').getCMSSWString().split(',')
#    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#    process.source.lumisToProcess.extend(myLumis)



if reReco == False:
    process.raw2digi_step = cms.Path(process.RawToDigi)
    process.L1Reco_step = cms.Path(process.L1Reco)

if hiReco:
    if reReco == False:
        process.reconstruction_step = cms.Path(process.reconstructionHeavyIons_withPF)
    else:
        process.rechits = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)
        process.hiTrackReco = cms.Sequence(process.rechits * process.heavyIonTracking)

        process.load("RecoHI.HiEgammaAlgos.HiElectronSequence_cff")
        process.load("RecoHI.HiEgammaAlgos.HiEgamma_cff")
        process.reconstruction_step =  cms.Path(
            process.hiTrackReco
            *process.muonRecoPbPb
            *process.hiPhotonSequence
            *process.electronGsfTrackingHi
            *process.HiParticleFlowLocalReco
            *process.HiParticleFlowReco
            )
        process.HiParticleFlowReco.remove(process.electronsWithPresel)
        process.HiParticleFlowReco.remove(process.electronsCiCLoose)
        
    #iterative tracking
    process.load("RecoHI.HiTracking.hiIterTracking_cff")
    process.heavyIonTracking *= process.hiIterTracking 

    # redo muons, seeding with global iterative tracks
    process.globalMuons.TrackerCollectionLabel = "hiGeneralTracks"
    process.muons.TrackExtractorPSet.inputTrackCollection = "hiGeneralTracks"
    process.muons.inputCollectionLabels = ["hiGeneralTracks", "globalMuons", "standAloneMuons:UpdatedAtVtx", "tevMuons:firstHit", "tevMuons:picky", "tevMuons:dyt"]

    # PF configuration
    process.particleFlowTmp.postMuonCleaning = False
    process.particleFlowBlock.RecMuons = 'muons'
    process.particleFlowClusterPS.thresh_Pt_Seed_Endcap = cms.double(99999.)
    process.pfTrack.UseQuality = True 
    process.pfTrack.TrackQuality = cms.string('highPurity')    
    process.pfTrack.TkColList = cms.VInputTag("hiGeneralTracks")
    process.pfTrack.GsfTracksInEvents = cms.bool(False)
    process.particleFlowTmp.usePFElectrons = False
    process.particleFlowTmp.useEGammaElectrons = False
    #process.particleFlowBlock.useEGPhotons = False
    process.pfTrackElec.applyGsfTrackCleaning = cms.bool(False)
    
else:
    process.reconstruction_step = cms.Path(process.reconstruction)


# clone one jet reco from the other
if hiReco:
    process.load("CmsHi.JetAnalysis.ExtraJetReco_cff")
    process.akPu3PFJets.src = 'particleTowerProducer'
    process.ak5PFJets = process.akPu3PFJets.clone(
        jetType = 'PFJet',
        src = 'particleFlowTmp',
        rParam = 0.5
        )

else:

    process.ak5PFJets.doAreaFastjet = False
    process.ak5PFJets.jetPtMin = 10.
    # 'Pu' in the name, just means towered
    process.akPu3PFJets = process.ak5PFJets.clone(
        jetType = 'BasicJet',
        src = 'particleTowerProducer',
        rParam = 0.3
        )

if hiMode:
    process.akPu3PFJets.doPUOffsetCorr = True
else:
    process.akPu3PFJets.doPUOffsetCorr = False
    process.ak5PFJets.doPUOffsetCorr = False

#process.akPu3PFJets.jetPtMin =1.
#process.ak5PFJets.jetPtMin =1.

process.load('RecoHI.HiJetAlgos.ParticleTowerProducer_cff')


process.reco_extra_jet    = cms.Path(         
    process.particleTowerProducer*
    process.akPu3PFJets
    )

if hiMode == False:
    process.reco_extra_jet *= process.ak5PFJets


# add gen step for HI

if isMC:
    process.load('CmsHi.JetAnalysis.ExtraGenReco_cff')
    process.HiGenParticleAna = cms.EDAnalyzer("HiGenAnalyzer")
    process.HiGenParticleAna.src= cms.untracked.InputTag(genTag)    
    process.hiGenParticles.srcVector = cms.vstring(genTag)
    #print "excluding neutrinos"
    #process.hiGenParticlesForJets.ignoreParticleIDs += cms.vuint32( 12,14,16)
    process.higen_step          = cms.Path(     
        process.hiGenParticles * process.hiGenParticlesForJets * process.genPartons * process.hiPartons * process.hiRecoGenJets #* process.HiGenParticleAna
        )
    
########  ADD EXTRA RECO WITH REGIT


if hiReco:
    if hiMode: svJetSel = 'pt > 30. && eta > -2. && eta < 2'
    else: svJetSel = 'pt > 10. && eta > -2. && eta < 2'

    print "Seeding around jets with, ",svJetSel

    from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
    # L2 (relative eta-conformity) Correction Services
    process.akPu3PFL2Relative = cms.ESSource(
        'LXXXCorrectionService',
        level     = cms.string('L2Relative'),
        algorithm = cms.string('AK3PF'),
        #era    = cms.string('HI_PFTowers_hiGoodTightTracks_D6T_413'),
        era    = cms.string('JEC_dijet'),
        section   = cms.string('')
        )
    # L3 (absolute) Correction Services
    process.akPu3PFL3Absolute = cms.ESSource(
        'LXXXCorrectionService',
        level     = cms.string('L3Absolute'),
        algorithm = cms.string('AK3PF'),
        #era    = cms.string('HI_PFTowers_hiGoodTightTracks_D6T_413'),
        era    = cms.string('JEC_dijet'),
        section   = cms.string('')
        )
    # L2L3 CORRECTION SERVICES
    process.akPu3PFL2L3 = cms.ESSource(
        'JetCorrectionServiceChain',
        correctors = cms.vstring('akPu3PFL2Relative','akPu3PFL3Absolute')
        )
    
    process.akPu3PFJetsL2L3   = cms.EDProducer('BasicJetCorrectionProducer',
                                               src         = cms.InputTag('akPu3PFJets'),
                                               correctors  = cms.vstring('akPu3PFL2L3')
                                               )
    
    
    
    process.ptEtaFilteredJets = cms.EDFilter("CandViewSelector",
                                             src = cms.InputTag("akPu3PFJetsL2L3"),
                                             cut = cms.string(svJetSel)
                                             )
    
    process.akPu3PFSelectedJets = cms.EDFilter("LargestPtCandViewSelector",
                                               src = cms.InputTag("ptEtaFilteredJets"),
                                               maxNumber = cms.uint32(3)
                                               )
        
    process.load("RecoHI.HiTracking.hiRegitTracking_cff")
    
    process.hiRegitInitialStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")
    process.hiRegitLowPtTripletStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")
    process.hiRegitPixelPairStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")
    process.hiRegitDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")
    process.hiRegitMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")
    process.hiRegitMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("akPu3PFSelectedJets")


    if hiMode == False:  # open up region for pp
        process.hiRegitInitialStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 2.
        process.hiRegitLowPtTripletStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 2.
        process.hiRegitPixelPairStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 2.
        process.hiRegitDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 15.
        process.hiRegitMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.originHalfLength = 10.
        process.hiRegitMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.originHalfLength = 10.

    process.load("RecoHI.HiTracking.MergeRegit_cff")
    
    process.regGlobalMuons = process.globalMuons.clone(
        TrackerCollectionLabel = "hiGeneralAndRegitTracks"
        )
    process.regGlbTrackQual = process.glbTrackQual.clone(
        InputCollection = "regGlobalMuons",
        InputLinksCollection = "regGlobalMuons"
        )
    process.regMuons = process.muons.clone()
    process.regMuons.TrackExtractorPSet.inputTrackCollection = "hiGeneralAndRegitTracks"
    process.regMuons.globalTrackQualityInputTag = "regGlbTrackQual"
    process.regMuons.inputCollectionLabels = cms.VInputTag("hiGeneralAndRegitTracks", "regGlobalMuons", "standAloneMuons:UpdatedAtVtx", "tevMuons:firstHit", "tevMuons:picky",
                                                           "tevMuons:dyt")
    
    
    process.regMuonReco = cms.Sequence(
        process.regGlobalMuons*
        process.regGlbTrackQual*
        process.regMuons
        )
        
    process.pfRegTrack = process.pfTrack.clone(TkColList = cms.VInputTag("hiGeneralAndRegitTracks"))
    process.pfRegBlock = process.particleFlowBlock.clone(RecTracks = "pfRegTrack",useIterTracking=True)
    process.regParticleFlow = process.particleFlowTmp.clone( blocks = "pfRegBlock")
    
    process.pfRegTrack.TrackQuality = "highPurity"
    
    process.pfRegTrack.GsfTracksInEvents = False
    process.regParticleFlow.usePFElectrons = False
    process.pfRegTrack.MuColl = "regMuons"
    process.pfRegBlock.RecMuons = "regMuons"
    process.regParticleFlow.muons = "regMuons"
    
    process.hiRegPF =  cms.Sequence(
        process.pfRegTrack
        *process.pfRegBlock
        *process.regParticleFlow
    )
    
    process.load("RecoHI.HiTracking.HICaloCompatibleTracks_cff")
    process.hiGeneralAndRegitCaloMatchedTracks = process.hiCaloCompatibleTracks.clone(    
        srcTracks = 'hiGeneralAndRegitTracks',
        srcPFCands = 'regParticleFlow'
        )

    svTrackSel = 'quality("highPurity") && pt > 1'
    print "Cutting on highPurity tracks"
    #if hiMode: svTrackSel = 'quality("highPurity") && pt > 1'
    #else: svTrackSel = 'quality("loose") && pt > 1'

    process.hiSecondaryVertexSelectedTracks = cms.EDFilter("TrackSelector",
                                                           src = cms.InputTag("hiGeneralAndRegitCaloMatchedTracks"),
                                                           cut = cms.string(
        svTrackSel
        )      
                                                           )

    # Need to re-reco PV with the same track collection as used for the SVs
    process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")
    process.offlinePrimaryVertices.TrackLabel = "hiSecondaryVertexSelectedTracks"


    process.regionalTracking = cms.Path(
        process.akPu3PFJetsL2L3 *
        process.ptEtaFilteredJets *
        process.akPu3PFSelectedJets *
        process.hiRegitTracking *
        process.hiGeneralAndRegitTracks *
        process.regMuonReco *
        process.hiRegPF *
        process.hiGeneralAndRegitCaloMatchedTracks *
        process.hiSecondaryVertexSelectedTracks *
        process.offlinePrimaryVertices
        )
    
    if redoPFJets:
        print "Redoing PF jets with regit PF "
        process.pfTowersFromRegit = process.particleTowerProducer.clone(src="regParticleFlow")
        process.akPu3regPFJets = process.akPu3PFJets.clone(src="pfTowersFromRegit")
        process.regionalTracking *= process.pfTowersFromRegit*process.akPu3regPFJets

# here comes pat
process.load('CmsHi.JetAnalysis.PatAna_cff')

# wow pat reconfigures the reco, how annoying
process.interestingEcalDetIdEB.basicClustersLabel = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters")
process.interestingEcalDetIdEE.basicClustersLabel = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters")

# Jet energy correction
process.jec = cms.ESSource("PoolDBESSource",
	DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
	timetype = cms.string('runnumber'),
	toGet = cms.VPSet(
		cms.PSet(record = cms.string("JetCorrectionsRecord"),
                         #tag = cms.string("JetCorrectorParametersCollection_HI_Calo_hiGoodTightTracks_D6T_413_IC5Calo"),
                         tag = cms.string('JetCorrectorParametersCollection_HI_Calo_hiGoodTightTracks_D6T_413_IC5Calo'),
                         label = cms.untracked.string("IC5Calo")),
                
		cms.PSet(record = cms.string("JetCorrectionsRecord"),
                         #tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_hiGoodTightTracks_D6T_413_AK3PF'),
                         tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_hiGoodTightTracks_PythiaZ2_442p5_AK3PF'),
                         label = cms.untracked.string("AK3PF")),
                
                ),
                           #connect = cms.string("sqlite_file:JEC_HI_PFTower_413patch2_2011_v3.db"),
                           connect = cms.string("sqlite_file:JEC_HI2760_CMSSW501_2012.db"),
                           )
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.akPu3PFcorr.payload = cms.string('AK3PF')

process.akPu3PFJetTracksAssociatorAtVertex.tracks = cms.InputTag(svTracks)
print "Restricting cone size to 0.3 "
process.akPu3PFJetTracksAssociatorAtVertex.coneSize=0.3
process.akPu3PFSecondaryVertexTagInfos.vertexCuts.maxDeltaRToJetAxis=0.3
process.akPu3PFImpactParameterTagInfos.primaryVertex = pvProducer

process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityComputer_cfi")
# these do nothing, jet probability is already positive
#process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityComputer_cfi")
#process.load("RecoBTag.ImpactParameterLearning.positiveOnlyJetBProbabilityComputer_cfi")  # Not in 44X
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D2ndComputer_cfi")
process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D3rdComputer_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeES_cfi")
process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexPositiveES_cfi")
process.load("RecoBTag.SoftLepton.negativeSoftLeptonByPtES_cfi")
process.load("RecoBTag.SoftLepton.positiveSoftLeptonByPtES_cfi")

# selection already done for SV
if hiReco == False:
    process.akPu3PFSecondaryVertexTagInfos.trackSelection.qualityClass = 'highPurity'
    process.akPu3PFSecondaryVertexNegativeTagInfos.trackSelection.qualityClass = 'highPurity'

process.akPu3PFpatJets.addAssociatedTracks = True
process.akPu3PFpatJets.addTagInfos = True
process.akPu3PFpatJets.addBTagInfo         = True
process.akPu3PFpatJets.addDiscriminators   = True
process.akPu3PFpatJets.getJetMCFlavour     = True


process.akPu3PFpatJets.tagInfoSources = cms.VInputTag(
    cms.InputTag("akPu3PFImpactParameterTagInfos"),
    cms.InputTag("akPu3PFSecondaryVertexTagInfos"),
    cms.InputTag("akPu3PFSoftMuonTagInfos"),
    )

# overloading akPu5 with ak5
process.akPu5PFcorr.src = "ak5PFJets"
process.akPu5PFJetTracksAssociatorAtVertex.jets   = "ak5PFJets"
process.akPu5PFSoftMuonTagInfos.jets           = "ak5PFJets"
process.akPu5PFPatJetPartonAssociation.jets    = "ak5PFJets"
process.akPu5PFmatch.src      = "ak5PFJets"
process.akPu5PFparton.src      = "ak5PFJets"
process.akPu5PFpatJets.jetSource            = "ak5PFJets"
                                                                     
process.akPu5PFcorr.payload = cms.string('AK5PF')

process.akPu5PFJetTracksAssociatorAtVertex.tracks = cms.InputTag(svTracks)
process.akPu5PFImpactParameterTagInfos.primaryVertex = pvProducer


if hiReco == False: process.akPu5PFSecondaryVertexTagInfos.trackSelection.qualityClass = 'highPurity'
    

process.akPu5PFpatJets.addAssociatedTracks = True
process.akPu5PFpatJets.addTagInfos = True
process.akPu5PFpatJets.addBTagInfo         = True
process.akPu5PFpatJets.addDiscriminators   = True
process.akPu5PFpatJets.getJetMCFlavour     = True


process.akPu5PFpatJets.tagInfoSources = cms.VInputTag(
    cms.InputTag("akPu5PFImpactParameterTagInfos"),
    cms.InputTag("akPu5PFSecondaryVertexTagInfos"),
    cms.InputTag("akPu5PFSoftMuonTagInfos"),
    )

if redoPFJets:
    print "resetting pat sequence to run regit jets"
    process.akPu3PFpatJets.jetSource = cms.InputTag("akPu3regPFJets")
    process.akPu3PFcorr.src= cms.InputTag("akPu3regPFJets")
    process.akPu3PFmatch.src= cms.InputTag("akPu3regPFJets")
    process.akPu3PFparton.src= cms.InputTag("akPu3regPFJets")
    process.akPu3PFPatJetPartonAssociation.jets    = cms.InputTag("akPu3regPFJets")
    process.akPu3PFJetTracksAssociatorAtVertex.jets   = cms.InputTag("akPu3regPFJets")
    process.akPu3PFSoftMuonTagInfos.jets = cms.InputTag("akPu3regPFJets")

if hiMode:
    if isMC: process.pat_step          = cms.Path(process.akPu3PFpatSequence_withBtagging)
    else: process.pat_step          = cms.Path(process.akPu3PFpatSequence_withBtagging_data)    
else:
    if isMC: process.pat_step          = cms.Path(process.akPu3PFpatSequence_withBtagging * process.akPu5PFpatSequence_withBtagging)
    else: process.pat_step          = cms.Path(process.akPu3PFpatSequence_withBtagging_data * process.akPu5PFpatSequence_withBtagging_data)



#######################
#   Analyzers
########################


process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(
    #'bTagAnalyzers.root'
    ivars.output
    )
                                   )


# Jet stuff

process.load('CmsHi.JetAnalysis.JetAnalyzers_cff')
process.akPu3PFJetAnalyzer.saveBfragments = True
if hiReco: process.akPu3PFJetAnalyzer.pfCandidateLabel = 'regParticleFlow'  #in order to get regit muons
process.akPu3PFJetAnalyzer.eventInfoTag = genTag
process.akPu3PFJetAnalyzer.hltTrgResults = cms.untracked.string(trigResults)
process.akPu3PFJetAnalyzer.useVtx = cms.untracked.bool(True)
process.akPu3PFJetAnalyzer.vtxTag = pvProducer
process.akPu3PFJetAnalyzer.trackTag = svTracks
if hiMode: process.akPu3PFJetAnalyzer.useCentrality = cms.untracked.bool(True)
else: process.akPu3PFJetAnalyzer.useCentrality = cms.untracked.bool(False)

if hiMode: process.akPu3PFJetAnalyzer.hltTrgNames = cms.untracked.vstring('HLT_HIMinBiasBSC_OR_v1','HLT_HIMinBiasHF_v1','HLT_HIMinBiasHf_OR_v1','HLT_HIMinBiasHfOrBSC_v1','HLT_HIL2Mu3_v1','HLT_HIL2Mu3_NHitQ_v1','HLT_HIL2Mu7_v1','HLT_HIL2Mu15_v1','HLT_HIJet55_v1','HLT_HIJet65_v1','HLT_HIJet80_v1','HLT_HIJet95_v1')
elif isMC: process.akPu3PFJetAnalyzer.hltTrgNames = cms.untracked.vstring('HLT_Jet15U_v3','HLT_Jet30U_v3','HLT_Jet50U_v3','HLT_Jet70U_v3','HLT_L2Mu7_v1','HLT_Mu0_v2','HLT_Mu3_v2','HLT_Mu5','HLT_Mu7','HLT_Mu9','HLT_L1_BptxXOR_BscMinBiasOR')
else:  process.akPu3PFJetAnalyzer.hltTrgNames = cms.untracked.vstring('HLT_HIMinBiasHfOrBSC_Core','HLT_Jet20_v1','HLT_Jet40_v1','HLT_Jet60_v1','HLT_L1SingleMuOpen_v1','HLT_Mu0_v3','HLT_Mu3_v3','HLT_Mu5_v3','HLT_L1BscMinBiasORBptxPlusANDMinus_v1')



if isMC:
    process.akPu3PFJetAnalyzer.isMC = True
else:
    process.akPu3PFJetAnalyzer.isMC = False
        
process.akPu3PFJetAnalyzer.doLifeTimeTagging = True
process.akPu3PFJetAnalyzer.jetTag = "akPu3PFpatJets"
process.akPu3PFJetAnalyzer.matchTag = "akPu5PFpatJets"


process.ak5PFJetAnalyzer = process.akPu3PFJetAnalyzer.clone(
    jetTag = "akPu5PFpatJets",
    genjetTag = "ak5HiGenJets",
    matchTag = "akPu3PFpatJets",
    rParam=0.5
    )


#  Track Analyzers

process.load("MitHig.PixelTrackletAnalyzer.trackAnalyzer_cff")

process.anaTrack.trackSrc = svTracks
process.anaTrack.vertexSrc = cms.vstring(pvProducer)
process.anaTrack.qualityString = "highPurity"

process.anaTrack.trackPtMin = 0
process.anaTrack.useQuality = False
process.anaTrack.doPFMatching = True
process.anaTrack.doSimTrack = hasSimInfo and isMC
if hiMode:process.anaTrack.useCentrality = True
else: process.anaTrack.useCentrality = False

process.load("SimTracker.TrackHistory.TrackClassifier_cff")

from SimTracker.TrackHistory.CategorySelectors_cff import *
process.hiTrackingCategorySelector = TrackingParticleCategorySelector(
    src = cms.InputTag('cutsTPForFak'),
    cut = cms.string("is('Bottom')")
)

process.hiTrackingCategorySelector.trackProducer = svTracks
process.hiTrackingCategorySelector.hepMC = cms.untracked.InputTag(genTag)
                            
process.bWeakDecaySelector = process.hiTrackingCategorySelector.clone(
    cut = cms.string("is('BWeakDecay')")
    )
process.cWeakDecaySelector = process.hiTrackingCategorySelector.clone(
    cut = cms.string("is('CWeakDecay')")
    )
process.v0DecaySelector = process.hiTrackingCategorySelector.clone(
    cut = cms.string("is('KsDecay') || is('LambdaDecay')")
    )
process.secondarySelector = process.hiTrackingCategorySelector.clone(
    cut = cms.string("is('SecondaryVertex')")
    )
process.tertiarySelector = process.hiTrackingCategorySelector.clone(
    cut = cms.string("is('TertiaryVertex')")
    )


process.bWeakDecayTracks = process.anaTrack.clone(
    doTrack = False,
    tpFakeSrc =  'bWeakDecaySelector',
    tpEffSrc = 'bWeakDecaySelector'
    )
process.cWeakDecayTracks = process.anaTrack.clone(
    doTrack = False,
    tpFakeSrc =  'cWeakDecaySelector',
    tpEffSrc = 'cWeakDecaySelector'
    )
process.v0DecayTracks = process.anaTrack.clone(
    doTrack = False,
    tpFakeSrc =  'v0DecaySelector',
    tpEffSrc = 'v0DecaySelector'
    )
process.secondaryTracks = process.anaTrack.clone(
    doTrack = False,
    tpFakeSrc =  'secondarySelector',
    tpEffSrc = 'secondarySelector'
    )
process.tertiaryTracks = process.anaTrack.clone(
    doTrack = False,
    tpFakeSrc =  'tertiarySelector',
    tpEffSrc = 'tertiarySelector'
    )

process.moreTrackAna = cms.Sequence(
    )


process.load("edwenger.HiTrkEffAnalyzer.hitrkEffAnalyzer_cff")

if hasSimInfo and isMC:                                                                                         
    process.trackAnalyzers = cms.Sequence(
        process.cutsTPForEff*
        process.cutsTPForFak*
        process.anaTrack*
        process.bWeakDecaySelector*process.bWeakDecayTracks*
        process.cWeakDecaySelector*process.cWeakDecayTracks*
        process.v0DecaySelector*process.v0DecayTracks*
        process.secondarySelector*process.secondaryTracks*
        process.tertiarySelector*process.tertiaryTracks
        )            
else:
    process.trackAnalyzers = cms.Sequence(
        process.anaTrack
        )            

# Muons
process.load("HiMuonAlgos.HLTMuTree.hltMuTree_cfi")
process.muonTree = process.hltMuTree.clone()
process.muonTree.doGen = cms.untracked.bool(True)
process.muonTree.vertices = pvProducer

process.ana_step          = cms.Path(         
    process.hiCentrality *
    process.akPu3PFJetAnalyzer *
    #process.trackAnalyzers*
    process.muonTree
    )

if hiMode==False: process.ana_step *= process.ak5PFJetAnalyzer

#IP calibration
process.load("RecoBTag.ImpactParameterLearning.ImpactParameterCalibration_cfi")
process.ipCalib.Jets                     = cms.InputTag('akPu3PFSelectedJets')
process.ipCalib.jetTagsColl              = cms.string("akPu3PFJetProbabilityBJetTags")
process.ipCalib.tagInfoSrc               = cms.InputTag("akPu3PFImpactParameterTagInfos")
process.ipCalib.jetPModuleName           = cms.string('akPu3PFJetProbabilityBJetTags')
process.ipCalib.jobIndex = ivars.jobIndex

#process.akPu3PFjetsIPcalib = cms.Path(process.ipCalib)

# trigger requirment, gets added in front of all patch w/ the superFilter
process.load('CmsHi.JetAnalysis.EventSelection_cff')
process.primaryVertexFilter.src = pvProducer
    
process.pvSel = cms.Path(process.primaryVertexFilter)
process.hbheNoiseSel = cms.Path(process.HBHENoiseFilter)
process.spikeSel = cms.Path(process.hiEcalRecHitSpikeFilter)
process.collSell = cms.Path(process.collisionEventSelection)
#process.hcalTimingSel = cms.Path(process.hcalTimingFilter)

process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.icPu5CaloJetsL2L3 = process.ak5CaloJetsL2L3.clone(src = 'iterativeConePu5CaloJets', correctors = ['ak5CaloL2L3'])

process.btagMuonInJet = cms.EDFilter("BTagSkimLeptonJet",
                                     CaloJet = cms.InputTag("icPu5CaloJetsL2L3"),
                                     MinimumCaloJetPt = cms.double(20.0),
                                     MinimumPtRel = cms.double(0.0),
                                     LeptonType = cms.string('muon'),
                                     Lepton = cms.InputTag("muons"),
                                     MinimumNLeptonJet = cms.int32(1),
                                     MaximumDeltaR = cms.double(0.4),
                                     MaximumLeptonEta = cms.double(2.5),
                                     MinimumLeptonPt = cms.double(5),
                                     MaximumCaloJetEta = cms.double(2.2)
                                     )
process.btagSkimPath       = cms.Path( process.icPu5CaloJetsL2L3 * process.btagMuonInJet)

process.load('CmsHi.HiHLTAlgos.hltanalysis_cff')
process.hltAna = cms.Path(process.hltanalysis)
process.pAna = cms.EndPath(process.skimanalysis)


process.hltanalysis.hltresults = cms.InputTag("TriggerResults")
process.skimanalysis.hltresults = cms.InputTag("TriggerResults")

# Secondary vertex sim matching, only tested on RAW data
'''
process.load("SimTracker.VertexAssociation.VertexAssociatorByTracks_cfi")

process.VertexAssociatorByTracksESProducer.trackingParticleSelector.stableOnlyTP= cms.bool(False)
if hiReco: process.VertexAssociatorByTracksESProducer.trackQuality = 'loose'  # already filtered on quality bit
else: process.VertexAssociatorByTracksESProducer.trackQuality = 'highPurity'  
    
# not sure what playback is for, maybe when RAWDEBUG is not stored
#process.load("SimTracker.TrackHistory.Playback_cff")
process.load("SimTracker.TrackHistory.SecondaryVertexTagInfoProxy_cff")
process.load("SimTracker.TrackHistory.VertexClassifier_cff")
process.svTagInfoProxy.svTagInfoProducer = cms.untracked.InputTag("akPu3PFSecondaryVertexTagInfos")

process.svTagInfoValidationAnalyzer = cms.EDAnalyzer("recoBSVTagInfoValidationAnalyzer",
                                                     process.vertexClassifier,
                                                     svTagInfoProducer = cms.untracked.InputTag('akPu3PFSecondaryVertexTagInfos'),
                                                     )
process.svTagInfoValidationAnalyzer.trackProducer = cms.untracked.InputTag(svTracks)
process.svTagInfoValidationAnalyzer.enableSimToReco = cms.untracked.bool(True)
process.svTagInfoValidationAnalyzer.enableRecoToSim = cms.untracked.bool(True)
process.svTagInfoValidationAnalyzer.vertexProducer = cms.untracked.InputTag('svTagInfoProxy')

process.svTagInfoSelector = SecondaryVertexTagInfoCategorySelector(
                                                                       src = cms.InputTag('akPu3PFSecondaryVertexTagInfos'),
                                                                       pxy = cms.InputTag('svTagInfoProxy'),
                                                                       #cut = cms.string("is('BWeakDecay') && !is('CWeakDecay')")
                                                                       #cut = cms.string("is('Fake')")
                                                                       cut = cms.string("is('SecondaryVertex')")
                                                                       )

process.svTagInfoProxy2 = cms.EDProducer('SecondaryVertexTagInfoProxy',
                                         svTagInfoProducer = cms.untracked.InputTag('svTagInfoSelector')
                                         )

process.vertexHistoryAnalyzer = cms.EDAnalyzer('VertexHistoryAnalyzer',
                                               process.vertexClassifier
                                               )
process.vertexHistoryAnalyzer.vertexProducer = 'svTagInfoProxy2'

process.svTagInfoSelector.trackProducer = cms.untracked.InputTag(svTracks)
process.vertexHistoryAnalyzer.trackProducer = cms.untracked.InputTag(svTracks)

process.svAna = cms.Path(#process.playback *
    process.svTagInfoProxy *
    process.svTagInfoValidationAnalyzer*
    process.svTagInfoSelector * process.svTagInfoProxy2 * process.vertexHistoryAnalyzer
    )
'''

if hltFilter:
    process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
    process.hltJetHI.HLTPaths = [hltFilter]
    process.hltJetHI.TriggerResultsTag = cms.InputTag(trigResults)
    process.superFilterSequence = cms.Sequence(process.hltJetHI)
    process.superFilterPath = cms.Path(process.superFilterSequence) 
    process.skimanalysis.superFilters = cms.vstring("superFilterPath")
    for path in process.paths:
        getattr(process,path)._seq = process.superFilterSequence*getattr(process,path)._seq

    process.hltAna.remove(process.hltJetHI)
    process.superFilterPath.remove(process.hltJetHI)
