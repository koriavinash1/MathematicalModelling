
//<![CDATA[

// a few things don't have var in front of them - they update already existing variables the game needs
lanesSide = 0;
patchesAhead = 1;
patchesBehind = 0;
trainIterations = 10000;

// the number of other autonomous vehicles controlled by your network
otherAgents = 0; // max of 10

var num_inputs = (lanesSide * 2 + 1) * (patchesAhead + patchesBehind);
var num_actions = 5;
var temporal_window = 3;
var network_size = num_inputs * temporal_window + num_actions * temporal_window + num_inputs;

var layer_defs = [];
    layer_defs.push({
    type: 'input',
    out_sx: 1,
    out_sy: 1,
    out_depth: network_size
});
layer_defs.push({
    type: 'fc',
    num_neurons: 1,
    activation: 'relu'
});
layer_defs.push({
    type: 'regression',
    num_neurons: num_actions
});

var tdtrainer_options = {
    learning_rate: 0.001,
    momentum: 0.0,
    batch_size: 64,
    l2_decay: 0.01
};

var opt = {};
opt.temporal_window = temporal_window;
opt.experience_size = 3000;
opt.start_learn_threshold = 500;
opt.gamma = 0.7;
opt.learning_steps_total = 10000;
opt.learning_steps_burnin = 1000;
opt.epsilon_min = 0.0;
opt.epsilon_test_time = 0.0;
opt.layer_defs = layer_defs;
opt.tdtrainer_options = tdtrainer_options;

brain = new deepqlearn.Brain(num_inputs, num_actions, opt);

learn = function (state, lastReward) {
brain.backward(lastReward);
var action = brain.forward(state);

draw_net();
draw_stats();

return action;
}

//]]>

/*###########*/
if (brain) {
brain.value_net.fromJSON({"layers":[{"out_depth":19,"out_sx":1,"out_sy":1,"layer_type":"input"},{"out_depth":1,"out_sx":1,"out_sy":1,"layer_type":"fc","num_inputs":19,"l1_decay_mul":0,"l2_decay_mul":1,"filters":[{"sx":1,"sy":1,"depth":19,"w":{"0":0.31567045383223724,"1":0.5541162400242423,"2":-0.06618880420925391,"3":0.13553373080567355,"4":0.1739899650526784,"5":-0.3134474434858062,"6":-0.05438945647171427,"7":0.374049860158203,"8":-0.06605495982402472,"9":-0.5936614419864751,"10":-0.036933857470188965,"11":0.551631982258335,"12":-0.4375874302350659,"13":0.33520608688221026,"14":-0.33579435090628434,"15":-0.04222879099074362,"16":0.33899457670536026,"17":-0.33542618166715654,"18":-0.35270539621434027}}],"biases":{"sx":1,"sy":1,"depth":1,"w":{"0":0.1}}},{"out_depth":1,"out_sx":1,"out_sy":1,"layer_type":"relu"},{"out_depth":5,"out_sx":1,"out_sy":1,"layer_type":"fc","num_inputs":1,"l1_decay_mul":0,"l2_decay_mul":1,"filters":[{"sx":1,"sy":1,"depth":1,"w":{"0":-0.9823303723897566}},{"sx":1,"sy":1,"depth":1,"w":{"0":2.051719508339062}},{"sx":1,"sy":1,"depth":1,"w":{"0":-0.07473837669606338}},{"sx":1,"sy":1,"depth":1,"w":{"0":1.4644305338905306}},{"sx":1,"sy":1,"depth":1,"w":{"0":-0.8467877768302512}}],"biases":{"sx":1,"sy":1,"depth":5,"w":{"0":0,"1":0,"2":0,"3":0,"4":0}}},{"out_depth":5,"out_sx":1,"out_sy":1,"layer_type":"regression","num_inputs":5}]});
}