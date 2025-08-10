# functions that are not used inside of the package but are managed by extensions are defined here.


export
    ####Makie
    record_video,

    ####SIMD
    #simd_overlaptest,
    simd_neighbor_traverse,
    simdbuild_traverse_bvh,
    #testsimdbuild_traverse_bvh,
    simdtwocluster_proximitytest!,
    simdneighbor_traverse,
    simdTreeData,
    ####experimental SIMD methods
    #simd_overlaptest,
    exptsimd_neighbor_traverse,
    exptsimdbuild_traverse_bvh,
    #testsimdbuild_traverse_bvh,
    exptsimdtwocluster_proximitytest!,

    ####GPU
    #gpuSpheresBVHSpecs,
    #gpuGridKey,
    gpubvh_neighborlist,
    gpubvh_neighborlist!,
    #gpuTreeData!,
    gpubounding_volume_hierarchy!,
    gpuneighbor_traverse!,
    gpubuild_traverse_bvh,
    PointPrimitive,
    gpuproximity_test!,
    gpuoverlap_test


function record_video end

function simd_overlaptest end

function simd_neighbor_traverse end
function simdneighbor_traverse end
function simdTreeData end

function simdbuild_traverse_bvh end


function simdtwocluster_proximitytest! end

function exptsimd_neighbor_traverse end
function exptsimdbuild_traverse_bvh end



function gpuSpheresBVHSpecs end
function gpuGridKey end

#function gpuTreeData! end
function gpubounding_volume_hierarchy! end
function gpuneighbor_traverse! end

function PointPrimitive end
function gpuproximity_test! end
function gpuoverlap_test end

function gpubvh_neighborlist end
function gpubvh_neighborlist! end