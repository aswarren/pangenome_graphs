package pangenome;

import java.util.ArrayList;
import java.util.List;
import org.gephi.graph.api.HierarchicalGraph;
import org.gephi.layout.plugin.AbstractLayout;
import org.gephi.layout.spi.Layout;
import org.gephi.layout.spi.LayoutBuilder;
import org.gephi.layout.spi.LayoutProperty;
import org.gephi.layout.plugin.force.yifanHu.YifanHuLayout;
import org.gephi.layout.plugin.force.yifanHu.YifanHuProportional;
import org.gephi.layout.plugin.multilevel.MaximalMatchingCoarsening;
import org.gephi.layout.plugin.random.RandomLayout;
import org.openide.util.NbBundle;

/**
 *
 * @author Helder Suzuki <heldersuzuki@gephi.org>
 */
public class MultiSpecial extends AbstractLayout implements Layout {

    private HierarchicalGraph graph;
    private int level;
    private YifanHuLayout layout;
    private YifanHuProportional yifanHu;
    private MaximalMatchingCoarsening coarseningStrategy;
    private int minSize;
    private double minCoarseningRate;
    private float stepRatio;
    private float optimalDistance;
    private int quadTreeMaxLevel;
    private float barnesHutTheta;

    //Security
    private int initedView;

    public MultiSpecial(LayoutBuilder layoutBuilder,
            MaximalMatchingCoarsening coarsening) {
        super(layoutBuilder);
        this.coarseningStrategy = coarsening;
        //     this.yifanHu = new YifanHu();
        this.yifanHu = new YifanHuProportional();
    }

    public void initAlgo() {
        graph = graphModel.getHierarchicalGraphVisible();
        initedView = graph.getView().getViewId();
        setConverged(false);
        level = 0;

        while (true) {
            int graphSize = graph.getTopNodes().toArray().length;
            coarseningStrategy.coarsen(graph);
            level++;
            int newGraphSize = graph.getTopNodes().toArray().length;
            if (newGraphSize < getMinSize() || newGraphSize > graphSize * getMinCoarseningRate()) {
                break;
            }
        }

        Layout random = new RandomLayout(null, 1000);
        random.setGraphModel(graphModel);
        random.initAlgo();
        random.goAlgo();

        initYifanHu();
    }

    void initYifanHu() {
        layout = yifanHu.buildLayout();
        layout.setGraphModel(graphModel);
        layout.resetPropertiesValues();
        layout.setAdaptiveCooling(false);
        layout.setStepRatio(stepRatio);
        layout.setOptimalDistance(optimalDistance);
        layout.setBarnesHutTheta(barnesHutTheta);
        layout.setQuadTreeMaxLevel(quadTreeMaxLevel);
        layout.initAlgo();
    }

    public void goAlgo() {
        HierarchicalGraph newGraph = graphModel.getHierarchicalGraphVisible();
        if(newGraph.getView().getViewId()!=initedView) {
            setConverged(true);
            layout.endAlgo();
            endAlgo();
            return;
        }
        this.graph = newGraph;
        if (layout.canAlgo()) {
            layout.goAlgo();
        } else {
            layout.endAlgo();
            if (level > 0) {
                coarseningStrategy.refine(graph);
                level--;

                initYifanHu();
            } else {
                setConverged(true);
                layout = null;
            }
        }
    }

    public void endAlgo() {
        while (level > 0) {
            coarseningStrategy.refine(graph);
            level--;
        }
    }

    public void resetPropertiesValues() {
        setMinSize(3);
        setMinCoarseningRate(0.75d);
        setStepRatio(0.97f);
        setOptimalDistance(100f);
        setQuadTreeMaxLevel(10);
        setBarnesHutTheta(1.2f);
    }

    public LayoutProperty[] getProperties() {
        List<LayoutProperty> properties = new ArrayList<LayoutProperty>();
        final String MULTILEVEL_CATEGORY = "Multi-level";
        final String YIFANHU_CATEGORY = "Yifan Hu's properties";
        final String BARNESHUT_CATEGORY = "Barnes-Hut's properties";

        try {
            properties.add(LayoutProperty.createProperty(
                    this, Integer.class, 
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.minimumLevelSize.name"),
                    MULTILEVEL_CATEGORY,
                    "YifanHuMultiLevel.minimumLevelSize.name",
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.minimumLevelSize.desc"),
                    "getMinSize", "setMinSize"));
            properties.add(LayoutProperty.createProperty(
                    this, Double.class, 
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.minimumCoarseningRate.name"),
                    MULTILEVEL_CATEGORY,
                    "YifanHuMultiLevel.minimumCoarseningRate.name",
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.minimumCoarseningRate.desc"),
                    "getMinCoarseningRate", "setMinCoarseningRate"));

            properties.add(LayoutProperty.createProperty(
                    this, Float.class, 
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.stepRatio.name"),
                    YIFANHU_CATEGORY,
                    "YifanHuMultiLevel.stepRatio.name",
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.stepRatio.desc"),
                    "getStepRatio", "setStepRatio"));
            properties.add(LayoutProperty.createProperty(
                    this, Float.class, 
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.optimalDistance.name"),
                    YIFANHU_CATEGORY,
                    "YifanHuMultiLevel.optimalDistance.name",
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.optimalDistance.desc"),
                    "getOptimalDistance", "setOptimalDistance"));

            properties.add(LayoutProperty.createProperty(
                    this, Integer.class, 
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.quadtreeMaxLevel.name"),
                    BARNESHUT_CATEGORY,
                    "YifanHuMultiLevel.quadtreeMaxLevel.name",
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.quadtreeMaxLevel.desc"),
                    "getQuadTreeMaxLevel", "setQuadTreeMaxLevel"));
            properties.add(LayoutProperty.createProperty(
                    this, Float.class, 
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.theta.name"),
                    BARNESHUT_CATEGORY,
                    "YifanHuMultiLevel.theta.name",
                    NbBundle.getMessage(getClass(), "YifanHuMultiLevel.theta.desc"),
                    "getBarnesHutTheta", "setBarnesHutTheta"));
        } catch (Exception e) {
            e.printStackTrace();
        }
        return properties.toArray(new LayoutProperty[0]);
    }

    /**
     * @return the minSize
     */
    public Integer getMinSize() {
        return minSize;
    }

    /**
     * @param minSize the minSize to set
     */
    public void setMinSize(Integer minSize) {
        this.minSize = minSize;
    }

    /**
     * @return the minCoarseningRate
     */
    public Double getMinCoarseningRate() {
        return minCoarseningRate;
    }

    /**
     * @param minCoarseningRate the minCoarseningRate to set
     */
    public void setMinCoarseningRate(Double minCoarseningRate) {
        this.minCoarseningRate = minCoarseningRate;
    }

    /**
     * @return the stepRatio
     */
    public Float getStepRatio() {
        return stepRatio;
    }

    /**
     * @param stepRatio the stepRatio to set
     */
    public void setStepRatio(Float stepRatio) {
        this.stepRatio = stepRatio;
    }

    /**
     * @return the optimalDistance
     */
    public Float getOptimalDistance() {
        return optimalDistance;
    }

    /**
     * @param optimalDistance the optimalDistance to set
     */
    public void setOptimalDistance(Float optimalDistance) {
        this.optimalDistance = optimalDistance;
    }

    /**
     * @return the quadTreeMaxLevel
     */
    public Integer getQuadTreeMaxLevel() {
        return quadTreeMaxLevel;
    }

    /**
     * @param quadTreeMaxLevel the quadTreeMaxLevel to set
     */
    public void setQuadTreeMaxLevel(Integer quadTreeMaxLevel) {
        this.quadTreeMaxLevel = quadTreeMaxLevel;
    }

    /**
     * @return the barnesHutTheta
     */
    public Float getBarnesHutTheta() {
        return barnesHutTheta;
    }

    /**
     * @param barnesHutTheta the barnesHutTheta to set
     */
    public void setBarnesHutTheta(Float barnesHutTheta) {
        this.barnesHutTheta = barnesHutTheta;
    }

    public interface CoarseningStrategy {

        public void coarsen(HierarchicalGraph graph);

        public void refine(HierarchicalGraph graph);
    }
}

