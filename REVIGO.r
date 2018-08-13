

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0051291","protein heterooligomerization", 0.003,-3.606,-4.634, 3.114,-4.4510,0.816,0.000),
c("GO:0052697","xenobiotic glucuronidation", 0.000, 4.105, 3.928, 1.398,-12.8268,0.671,0.000),
c("GO:0009812","flavonoid metabolic process", 0.005,-3.833, 0.511, 3.375,-9.2097,0.933,0.021),
c("GO:0021549","cerebellum development", 0.005,-2.408,-2.040, 3.437,-3.8827,0.898,0.023),
c("GO:0009813","flavonoid biosynthetic process", 0.004,-3.218, 3.607, 3.311,-8.1543,0.878,0.025),
c("GO:0032776","DNA methylation on cytosine", 0.035, 4.719,-5.050, 4.238,-6.0540,0.613,0.028),
c("GO:0042440","pigment metabolic process", 0.490,-2.244, 5.551, 5.387,-5.2321,0.906,0.057),
c("GO:0009698","phenylpropanoid metabolic process", 0.006,-4.974, 2.121, 3.480,-4.3170,0.869,0.084),
c("GO:0006069","ethanol oxidation", 0.015, 0.165, 6.407, 3.886,-3.1701,0.879,0.103),
c("GO:0002230","positive regulation of defense response to virus by host", 0.001, 4.856,-2.929, 2.593,-3.6556,0.768,0.168),
c("GO:0010225","response to UV-C", 0.001, 0.219, 0.646, 2.494,-3.1261,0.917,0.187),
c("GO:0007213","G-protein coupled acetylcholine receptor signaling pathway", 0.003, 4.243,-0.651, 3.132,-3.0057,0.753,0.207),
c("GO:0046602","regulation of mitotic centrosome separation", 0.000, 1.539,-5.062, 1.839,-3.0320,0.735,0.251),
c("GO:0006109","regulation of carbohydrate metabolic process", 0.056, 6.546, 1.624, 4.441,-3.5171,0.688,0.266),
c("GO:0042573","retinoic acid metabolic process", 0.001, 5.778, 3.448, 2.669,-4.0039,0.723,0.282),
c("GO:0032787","monocarboxylic acid metabolic process", 2.312, 2.201, 6.174, 6.061,-3.1713,0.794,0.285),
c("GO:0040029","regulation of gene expression, epigenetic", 0.225, 6.552,-4.069, 5.049,-4.1029,0.714,0.306),
c("GO:0006063","uronic acid metabolic process", 0.027, 4.681, 5.494, 4.134,-9.7011,0.736,0.383),
c("GO:0006304","DNA modification", 0.364, 0.714,-7.600, 5.257,-3.8894,0.800,0.416),
c("GO:0042501","serine phosphorylation of STAT protein", 0.001, 3.790,-3.607, 2.577,-3.6021,0.696,0.420),
c("GO:0045922","negative regulation of fatty acid metabolic process", 0.018, 6.734, 0.684, 3.944,-9.7011,0.524,0.437),
c("GO:0006335","DNA replication-dependent nucleosome assembly", 0.001,-2.013,-5.826, 2.571,-3.7375,0.720,0.461),
c("GO:0045926","negative regulation of growth", 0.014, 8.202,-1.826, 3.835,-3.0343,0.693,0.494),
c("GO:0010677","negative regulation of cellular carbohydrate metabolic process", 0.002, 7.703, 0.484, 2.923,-5.3665,0.592,0.499),
c("GO:1902037","negative regulation of hematopoietic stem cell differentiation", 0.000, 7.229,-3.003, 1.491,-3.0009,0.687,0.502),
c("GO:0005996","monosaccharide metabolic process", 1.214, 4.237, 6.224, 5.781,-3.0429,0.773,0.508),
c("GO:0050691","regulation of defense response to virus by host", 0.001, 4.575,-2.438, 2.821,-3.2984,0.763,0.512),
c("GO:0035235","ionotropic glutamate receptor signaling pathway", 0.011, 4.977,-0.920, 3.751,-3.0057,0.739,0.515),
c("GO:0016458","gene silencing", 0.025, 7.689,-2.418, 4.090,-5.7447,0.624,0.551),
c("GO:0043414","macromolecule methylation", 1.178, 1.776,-7.906, 5.768,-3.5114,0.795,0.562),
c("GO:0060968","regulation of gene silencing", 0.006, 8.071,-2.228, 3.490,-5.7328,0.628,0.567),
c("GO:0065003","macromolecular complex assembly", 0.676,-3.166,-5.068, 5.527,-3.5622,0.795,0.578),
c("GO:0009804","coumarin metabolic process", 0.001,-5.162, 1.376, 2.606,-4.3170,0.876,0.590),
c("GO:0051290","protein heterotetramerization", 0.001,-3.544,-4.518, 2.468,-3.9136,0.826,0.592),
c("GO:0010565","regulation of cellular ketone metabolic process", 0.043, 7.088,-0.192, 4.332,-3.3686,0.673,0.667),
c("GO:0006805","xenobiotic metabolic process", 0.067, 1.256, 1.636, 4.526,-5.6925,0.855,0.668),
c("GO:0045912","negative regulation of carbohydrate metabolic process", 0.002, 7.650, 0.743, 3.027,-5.2083,0.604,0.689));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);



# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");
