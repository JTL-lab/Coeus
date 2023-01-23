export default {
  plot: {
    transitionDuration: 250,
    scaleFactor: 15,
    scaleGenes: true,
    fontFamily:
      'system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Ubuntu, "Helvetica Neue", Oxygen, Cantarell, sans-serif',
  },
  legend: {
    entryHeight: 18,
    fontSize: 14,
    onClickCircle: null,
    onClickText: null,
    show: true,
    /* Control gene group legend placement in dashboard*/
    marginLeft: 50,
    marginTop: 100,
  },
  /* Controls % identity bar */
  colourBar: {
    fontSize: 10,
    height: 14,
    show: true,
    width: 150,
    marginTop: 50,
    marginLeft: 2200,
  },
  /* Controls scale bar next to % identity bar */
  scaleBar: {
    colour: "black",
    fontSize: 10,
    height: 12,
    basePair: 2500,
    show: true,
    stroke: 1,
    marginTop: 50,
    marginLeft: 2000,
  },
  link: {
    show: true,
    asLine: false,
    straight: false,
    threshold: 0,
    strokeWidth: 0.5,
    groupColour: false,
    bestOnly: false,
    label: {
      show: false,
      fontSize: 10,
      background: true,
      position: 0.5,
    },
  },
  cluster: {
    nameFontSize: 12,
    lociFontSize: 10,
    hideLocusCoordinates: false,
    spacing: 30,
    alignLabels: true,
  },
  locus: {
    trackBar: {
      colour: "#000000",
      stroke: 1,
    },
    spacing: 50,
  },
  gene: {
    shape: {
      bodyHeight: 12,
      tipHeight: 5,
      tipLength: 12,
      onClick: null,
      stroke: "black",
      strokeWidth: 1,
    },
    label: {
      anchor: "start",
      fontSize: 10,
      rotation: 25,
      position: "top",
      spacing: 2,
      show: false,
      start: 0.5,
      name: "uid",
    },
  },
};
