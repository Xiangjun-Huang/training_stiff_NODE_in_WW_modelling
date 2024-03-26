from NODE_model import neuralODE

def model_initiation(noFeature, yMu, yStd, dyMu, dyStd, yMax, yMin, dyMax, dyMin, ytScale):

    model = neuralODE(noFeature, yMu, yStd, dyMu, dyStd, yMax, yMin, dyMax, dyMin, ytScale)

    return model