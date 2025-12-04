def parseCliqueBaitJson(json_file) {
    def jsonSlurper = new groovy.json.JsonSlurper()
    def parsed = jsonSlurper.parseText(json_file.text)
    def outList = []
    tax_name = json_file.toString().split("/")[-1].replace("_baits.json","")
    for (item in parsed) {
        def tt = tuple(tax_name + "_" + item.id.toString(), item.genomes )
        outList << tt
    }
    return outList
}