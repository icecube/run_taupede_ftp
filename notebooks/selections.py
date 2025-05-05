

def create_selections( file ):

    selections = {}

    selections["RecoETot>60TeV"] = ( file['variables']["RecoETot"] > 60e3 )

    selections["TrueTrack"]         = selections["RecoETot>60TeV"] & ( file['variables']["MCInteractionEventclass"] == 3.0 )
    selections["TrueDoubleCascade"] = selections["RecoETot>60TeV"] & ( file['variables']["MCInteractionEventclass"] == 2.0 )
    selections["TrueSingleCascade"] = selections["RecoETot>60TeV"] & ( file['variables']["MCInteractionEventclass"] == 1.0 )

    selections["Track"]         = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass"] == 3.0 )
    selections["DoubleCascade"] = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass"] == 2.0 )
    selections["SingleCascade"] = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass"] == 1.0 )

    # selections["TrackEvtGen"]         = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass_evtgen"] == 3.0 )
    # selections["DoubleCascadeEvtGen"] = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass_evtgen"] == 2.0 )
    # selections["SingleCascadeEvtGen"] = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass_evtgen"] == 1.0 )

    
    selections["DoubleCascade_TrueTrack"]         = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass"] == 2.0 ) & ( file['variables']["MCInteractionEventclass"] == 3.0 )
    selections["DoubleCascade_TrueDoubleCascade"] = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass"] == 2.0 ) & ( file['variables']["MCInteractionEventclass"] == 2.0 )
    selections["DoubleCascade_TrueSingleCascade"] = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass"] == 2.0 ) & ( file['variables']["MCInteractionEventclass"] == 1.0 )

    # selections["DoubleCascadeEvtGen_TrueTrack"]         = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass_evtgen"] == 2.0 ) & ( file['variables']["MCInteractionEventclass"] == 3.0 )
    # selections["DoubleCascadeEvtGen_TrueDoubleCascade"] = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass_evtgen"] == 2.0 ) & ( file['variables']["MCInteractionEventclass"] == 2.0 )
    # selections["DoubleCascadeEvtGen_TrueSingleCascade"] = selections["RecoETot>60TeV"] & ( file['variables']["FinalEventClass_evtgen"] == 2.0 ) & ( file['variables']["MCInteractionEventclass"] == 1.0 )

    selections["SingleCascade_TrueTrack"]         = ( selections["RecoETot>60TeV"] ) & ( file['variables']["FinalEventClass"] == 1.0 ) & ( file['variables']["MCInteractionEventclass"] == 3.0 )
    selections["SingleCascade_TrueDoubleCascade"] = ( selections["RecoETot>60TeV"] ) & ( file['variables']["FinalEventClass"] == 1.0 ) & ( file['variables']["MCInteractionEventclass"] == 2.0 )
    selections["SingleCascade_TrueSingleCascade"] = ( selections["RecoETot>60TeV"] ) & ( file['variables']["FinalEventClass"] == 1.0 ) & ( file['variables']["MCInteractionEventclass"] == 1.0 )


    return selections
