/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow validate_parameters {
    
    main:
        Pinguscript.ping_start(nextflow, workflow, params)
        parameters_passed = true
        print("Checking parameters")

        available = Runtime.runtime.availableProcessors()
        
        if (params.threads as int > available){
            println("Specified number of threads ${params.threads} is greater than the number of available Processors ${available}")
            println("Number of threads is reset to ${available}")
        }

        if (params.human | params.pombe | params.cerevisiae){

        }
        else {
            parameters_passed = false
            println("Error - Species Must be Set")
        }

        try {
            file(params.modbam, checkIfExists:true)
        }
        catch (Exception e) {
            parameters_passed = false
            println("Error - ModBam Doesn't Exist")
        }

        try {
            file(params.telo_stats, checkIfExists:true)
        }
        catch (Exception e) {
            parameters_passed = false
            println("Error - Telo Stats File Doesn't Exist")
        }

        try {
            file(params.reference, checkIfExists:true)
        }
        catch (Exception e) {
            parameters_passed = false
            println("Error - Reference Doesn't Exist")
        }

        if (params.cluster_results != "") {
            try {
                file(params.cluster_results, checkIfExists:true)
            }
            catch (Exception e) {
                parameters_passed = false
                println("Error - Cluster File Doesn't Exist")
            }
        }

        if (params.spike_in_reference != ""){
            try {
                file(params.modbam, checkIfExists:true)
            }
            catch (Exception e) {
                parameters_passed = false
                println("Error - ModBam Doesn't Exist")
            }
        }

    emit:
        passed = parameters_passed
}
