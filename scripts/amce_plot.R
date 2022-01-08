library("cjoint")

amce_plot <- function(path.data, path.plot) {
    d.conjoint <- read.csv(path.data)
    d.conjoint$TECHNOLOGIE <- factor(
        d.conjoint$TECHNOLOGIE,
        levels=c("Photovoltaikanlagen auf Dächern", "Photovoltaik-Freiflächenanlagen", "Windturbinen an Land")
    )
    d.conjoint$MENGE_DER_STROMIMPORTE <- factor(
        d.conjoint$MENGE_DER_STROMIMPORTE,
        levels=c("Keine - ihr Strom stammt aus regionalen Anlagen.", "Niedrig", "Mittel", "Hoch")
    )
    d.conjoint$FL_CHENBEDARF <- factor(
        d.conjoint$FL_CHENBEDARF,
        levels=c("Kein weiterer Flächenbedarf", "Sehr gering", "Geringe", "Mittel", "Hoch", "Sehr hoch")
    )
    d.conjoint$STROMPREISENTWICKLUNG_F_R_HAUSHALTE <- factor(
        d.conjoint$STROMPREISENTWICKLUNG_F_R_HAUSHALTE,
        levels=c("Heutiges Niveau", "Minimaler Anstieg", "Moderater Anstieg", "Starker Anstieg")
    )
    d.conjoint$AUSBAU_DER_BERTRAGUNGSNETZ <- factor(
        d.conjoint$AUSBAU_DER_BERTRAGUNGSNETZ,
        levels=c("Heutiges Niveau", "Leichter Rückgang", "Leichter Anstieg", "Moderater Anstieg", "Starker Anstieg")
    )
    d.conjoint$EIGENTUMSVERH_LTNISSE_DER_ERZEUGUNGSANLAGEN <- factor(
        d.conjoint$EIGENTUMSVERH_LTNISSE_DER_ERZEUGUNGSANLAGEN,
        levels=c("Private Versorgungsunternehmen", "Öffentliche Träger", "Lokale und regionale Gemeinschaften")
    )
    results <- amce(CHOICE_INDICATOR ~  TECHNOLOGIE + MENGE_DER_STROMIMPORTE + FL_CHENBEDARF
        + STROMPREISENTWICKLUNG_F_R_HAUSHALTE + AUSBAU_DER_BERTRAGUNGSNETZ
        + EIGENTUMSVERH_LTNISSE_DER_ERZEUGUNGSANLAGEN,
        data=d.conjoint,
        cluster=TRUE,
        respondent.id="RESPONDENT_ID")
    png(path.plot)
    plot(
        results,
        xlab="Change in Pr(Design preferred)",
        xlim=c(-.5,.5),
        breaks=c(-.4, 0, .4),
        labels=c("-.4","0",".4"),
        text.size=13
    )
    dev.off()

}

amce_plot(snakemake@input[["data"]], snakemake@output[[1]])
